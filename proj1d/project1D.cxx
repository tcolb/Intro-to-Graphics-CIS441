#include <iostream>
#include <array>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;

enum TriangleType {GoingUp, GoingDown, Arbitrary};
typedef std::array<double, 3> tvec;
const int WINDOW_HEIGHT = 1000;
const int WINDOW_WIDTH = 1000;
int TRIANGLE_NUM = 0;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

class Triangle
{
  public:
      tvec X;
      tvec Y;
      tvec Z;
      std::array<tvec, 3> colors;
      TriangleType type;
      int vecL;
      int vecR;
      int vecD;

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      vtkImageData *image;
      unsigned char   *buffer;
      int width, height;
      //double **zbuf;  
      std::array<std::array<double, WINDOW_WIDTH>, WINDOW_HEIGHT> zbuf;

      Screen () {
        for (int i =0; i < WINDOW_HEIGHT; i++) {
            this->zbuf[i].fill(-1);
        }
      }

      ~Screen() { }

      void zbufferPrint() { 
        for (int i =0; i < WINDOW_HEIGHT; i++) {
            for (int j = 0; j < WINDOW_WIDTH; j++) {
                cout << this->zbuf[i][j] << " ";
            }
            cout << endl;
        }
      }

      double convertFRGB(double a) {
        return ceil_441(a * 255);
    }
   

void setPixel(double x, double y, double z, tvec rgb) {
           if (y < 0 || y >= WINDOW_HEIGHT || x < 0 || x >= WINDOW_WIDTH)
               return;

           if (this->zbuf[y][x] > z)
              return;

           this->zbuf[y][x] = z;
           this->buffer = (unsigned char*) image->GetScalarPointer(x, y , z);
           this->buffer[0] = convertFRGB(rgb[0]);
           this->buffer[1] = convertFRGB(rgb[1]);
           this->buffer[2] = convertFRGB(rgb[2]);
      }

      bool zbufferSet(double y, double x, double c) {
          if (y < 0 || y >= WINDOW_HEIGHT || x < 0 || x >= WINDOW_WIDTH)
              return false; 
          if (this->zbuf[y][x] < c) {
             this->zbuf[y][x] = c;
             return true;
          }
          return false;
      }
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double calcT(double x, double a, double b) {
    return double((x-a)/(b-a));   
}

double lerp(double t, double fa, double fb) {
    return double(fa + t*(fb-fa));
}

void color_lerp(double t, tvec aColor, tvec bColor, tvec &out ) {
    out[0] = lerp(t, aColor[0], bColor[0]);
    out[1] = lerp(t, aColor[1], bColor[1]);
    out[2] = lerp(t, aColor[2], bColor[2]);
}

void setTriangleProps(Triangle* t) {

    std::array<std::array<double, 6>, 3> verts;
    verts[0] = { t->X[0], t->Y[0], t->Z[0], t->colors[0][0], t->colors[0][1], t->colors[0][2] };
    verts[1] = { t->X[1], t->Y[1], t->Z[1], t->colors[1][0], t->colors[1][1], t->colors[1][2] };
    verts[2] = { t->X[2], t->Y[2], t->Z[2], t->colors[2][0], t->colors[2][1], t->colors[2][2] };

    std::sort(verts.begin(), verts.end(),
    [](const std::array<double, 6> &a1, const std::array<double, 6> &a2) {
        return (a1[1] < a2[1]);    
    });

    if (verts[0][1] == verts[1][1]) {
        t->type = GoingUp;
        if (verts[0][0] < verts[1][0]) {
            t->vecL = 0;
            t->vecR = 1;
        } else {
            t->vecL = 1;
            t->vecR = 0;
        }
        t->vecD = 2;
    } else if (verts[1][1] == verts[2][1]) {
        t->type = GoingDown;
        if (verts[1][0] < verts[2][0]) {
            t->vecL = 1;
            t->vecR = 2;
        } else {
            t->vecL = 2;
            t->vecR = 1;
        }
        t->vecD = 0;
    } else {
        t->type = Arbitrary;
        t->vecL = 0;
        t->vecR = 2;
        t->vecD = 1;
    }

    t->X[0] = verts[0][0];
    t->X[1] = verts[1][0];
    t->X[2] = verts[2][0];
    t->Y[0] = verts[0][1];
    t->Y[1] = verts[1][1];
    t->Y[2] = verts[2][1];
    t->Z[0] = verts[0][2];
    t->Z[1] = verts[1][2];
    t->Z[2] = verts[2][2];
    t->colors[0][0] = verts[0][3];
    t->colors[0][1] = verts[0][4];
    t->colors[0][2] = verts[0][5];
    t->colors[1][0] = verts[1][3];
    t->colors[1][1] = verts[1][4];
    t->colors[1][2] = verts[1][5];
    t->colors[2][0] = verts[2][3];
    t->colors[2][1] = verts[2][4];
    t->colors[2][2] = verts[2][5];

}

void rasterizeTriangle(Triangle tri, Screen* screen) {
    int vecD = tri.vecD;
    int vecL = tri.vecL;
    int vecR = tri.vecR;

    double rowMin, rowMax;
    
    if (tri.type == GoingDown) {
        rowMin = ceil_441(tri.Y[vecD]);
        rowMax = floor_441(tri.Y[vecL]);
        //cout << "Down ";
    } else if (tri.type == GoingUp) {
        rowMin = ceil_441(tri.Y[vecL]);
        rowMax = floor_441(tri.Y[vecD]);
        //cout << "Up ";
    }
    /*
    cout << "triangle: " << std::setprecision(16) << TRIANGLE_NUM << endl
         << "det: (" << tri.X[vecD] << ", " << tri.Y[vecD] << ", " << tri.Z[vecD] << ") "
         << ", color = (" << tri.colors[vecD][0] << ", " << tri.colors[vecD][1] << ", " << tri.colors[vecD][2] << endl
         << "left: (" << tri.X[vecL] << ", " << tri.Y[vecL] << ", " << tri.Z[vecL] << ") "
         << ", color = (" << tri.colors[vecL][0] << ", " << tri.colors[vecL][1] << ", " << tri.colors[vecL][2] << endl
         << "right: (" << tri.X[vecR] << ", " << tri.Y[vecR] << ", " << tri.Z[vecR] << ") " 
         << ", color = (" << tri.colors[vecR][0] << ", " << tri.colors[vecR][1] << ", " << tri.colors[vecR][2] << endl
         << "lineMin = " << rowMin << ", lineMax = " << rowMax << endl << endl;
    */

    // scanline row
    for (double row = rowMin; row <= rowMax; row++) {
        double t = calcT(row, tri.Y[vecD], tri.Y[vecL]);
        
        // not ceil/floored as using for color lerp as well
        double lEnd = lerp(t, tri.X[vecD], tri.X[vecL]);
        double rEnd = lerp(t, tri.X[vecD], tri.X[vecR]);

        // color lerping
        tvec lColor, rColor;
        color_lerp(t, tri.colors[vecD], tri.colors[vecL], lColor);
        color_lerp(t, tri.colors[vecD], tri.colors[vecR], rColor);
        // z-buffer lerping
        double lZ = lerp(t, tri.Z[vecD], tri.Z[vecL]);
        double rZ = lerp(t, tri.Z[vecD], tri.Z[vecR]);

        // scanline col
        for (double col = ceil_441(lEnd); col <= floor_441(rEnd); col++) {
            double xT = calcT(col, lEnd, rEnd);
            double xZ = lerp(xT, lZ, rZ);
            //cout << "zbuf access: " << row << col << endl;
            tvec xColor; 
            color_lerp(xT, lColor, rColor, xColor);
            screen->setPixel(col, row, xZ, xColor);
        }
    }
}

void rasterizeArbitraryTriangle(Triangle tri, Screen* screen) {
    int vecD = tri.vecD;
    int vecL = tri.vecL;
    int vecR = tri.vecR;

    Triangle upHalf = Triangle();
    Triangle downHalf = Triangle();

    double longT = calcT(tri.Y[vecD], tri.Y[vecL], tri.Y[vecR]);
    double longEnd = lerp(longT, tri.X[vecL], tri.X[vecR]);

/*
    cout << "ARB triangle: " << std::setprecision(16) << TRIANGLE_NUM << endl
         << "det: (" << tri.X[vecD] << ", " << tri.Y[vecD] << ", " << tri.Z[vecD] << ")" << endl
         << "left: (" << tri.X[vecL] << ", " << tri.Y[vecL] << ", " << tri.Z[vecL] << ")" << endl
         << "right: (" <creen->setPixel(col, row, xZ, xColor); tri.X[vecR] << ", " << tri.Y[vecR] << ", " << tri.Z[vecR] << ")" << endl << endl; 
*/

    // fixing right / left verts
    int newL = 0;
    int newR = 1;
    if (tri.X[vecD] < longEnd) { newL = 1; newR = 0; }

    // lerp endpoint colors for new tris
    tvec longColor;
    color_lerp(longT, tri.colors[vecL], tri.colors[vecR], longColor);

    // lerp endpoint z values for new tris
    double longZ = lerp(longT, tri.Z[vecL], tri.Z[vecR]);

    // setting going up
    upHalf.X[0] = longEnd;
    upHalf.X[1] = tri.X[vecD];
    upHalf.X[2] = tri.X[vecR];
    upHalf.Y[0] = tri.Y[vecD];
    upHalf.Y[1] = tri.Y[vecD];
    upHalf.Y[2] = tri.Y[vecR];
    upHalf.Z[0] = longZ;
    upHalf.Z[1] = tri.Z[vecD];
    upHalf.Z[2] = tri.Z[vecR];
    upHalf.vecL = newL;
    upHalf.vecR = newR;
    upHalf.vecD = 2;
    upHalf.type = GoingUp;
    upHalf.colors[0] = longColor;
    upHalf.colors[1] = tri.colors[vecD];
    upHalf.colors[2] = tri.colors[vecR];

    // setting going down
    downHalf.X[0] = longEnd;
    downHalf.X[1] = tri.X[vecD];
    downHalf.X[2] = tri.X[vecL];
    downHalf.Y[0] = tri.Y[vecD];
    downHalf.Y[1] = tri.Y[vecD];
    downHalf.Y[2] = tri.Y[vecL];
    downHalf.Z[0] = longZ;
    downHalf.Z[1] = tri.Z[vecD];
    downHalf.Z[2] = tri.Z[vecL];
    downHalf.vecL = newL;
    downHalf.vecR = newR;
    downHalf.vecD = 2;
    downHalf.type = GoingDown;
    downHalf.colors[0] = longColor;
    downHalf.colors[1] = tri.colors[vecD];
    downHalf.colors[2] = tri.colors[vecL];

    rasterizeTriangle(upHalf, screen);
    rasterizeTriangle(downHalf, screen);
}

void scanline(std::vector<Triangle> triangles, Screen* screen) {
    for (auto & t : triangles) {
        setTriangleProps(&t);

        if (t.type == Arbitrary) { rasterizeArbitraryTriangle(t, screen); }
        else { rasterizeTriangle(t, screen); }

        TRIANGLE_NUM++;
        //if (TRIANGLE_NUM == 10)
          //return;
    }
}

int main()
{
    vtkImageData *image = NewImage(WINDOW_WIDTH, WINDOW_HEIGHT);
    unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = WINDOW_WIDTH*WINDOW_HEIGHT;
    for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.image = image;
    screen.buffer = buffer;
    screen.width = WINDOW_WIDTH;
    screen.height = WINDOW_HEIGHT;

    scanline(triangles, &screen);
    
    WriteImage(image, "allTriangles");
}
