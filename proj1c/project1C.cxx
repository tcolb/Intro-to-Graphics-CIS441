#include <iostream>
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
const int WINDOW_HEIGHT = 1344;
const int WINDOW_WIDTH = 1786;
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
      double         X[3];
      double         Y[3];
      unsigned char color[3];
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

      void setPixel(double x, double y, double z, unsigned char* rgb) {
          if (y < 0 || y >= WINDOW_HEIGHT || x < 0 || x >= WINDOW_WIDTH)
              return;
          this->buffer = (unsigned char*) image->GetScalarPointer(x, y , z);
          this->buffer[0] = rgb[0];
          this->buffer[1] = rgb[1];
          this->buffer[2] = rgb[2];
      }  
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
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
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }
    cerr << "Done reading" << endl;

    return tris;
}

double scanline_slope(double x1, double y1, double x2, double y2) {
    if (x1 == x2)
        return 0;
    return double((y2-y1)/(x2-x1));
}

double scanline_b(double y, double m, double x) {
    return double(y-m*x);
}

double scanline_end(double x, double y, double m, double b) {
    if (m == 0)
        return x;
    return double((y-b)/m);
}

void setTriangleProps(Triangle* t) {
    std::vector< std::pair<double, double> > sortVec;
    sortVec.push_back(std::make_pair(t->X[0], t->Y[0]));
    sortVec.push_back(std::make_pair(t->X[1], t->Y[1]));
    sortVec.push_back(std::make_pair(t->X[2], t->Y[2]));
    std::sort(sortVec.begin(), sortVec.end(),
    [](const std::pair<double,double> &p1, const std::pair<double,double> &p2) {
        return (p1.second < p2.second);    
    });

    if (sortVec[0].second == sortVec[1].second) {
        t->type = GoingUp;
        if (sortVec[0].first < sortVec[1].first) {
            t->vecL = 0;
            t->vecR = 1;
        } else {
            t->vecL = 1;
            t->vecR = 0;
        }
        t->vecD = 2;
    } else if (sortVec[1].second == sortVec[2].second) {
        t->type = GoingDown;
        if (sortVec[1].first < sortVec[2].first) {
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

    t->X[0] = sortVec[0].first;
    t->X[1] = sortVec[1].first;
    t->X[2] = sortVec[2].first;
    t->Y[0] = sortVec[0].second;
    t->Y[1] = sortVec[1].second;
    t->Y[2] = sortVec[2].second;

}

void rasterizeTriangle(Triangle triangle, Screen screen) {
    int vecD = triangle.vecD;
    int vecL = triangle.vecL;
    int vecR = triangle.vecR;

    double lineMin, lineMax;
    
    if (triangle.type == GoingDown) {
        lineMin = ceil_441(triangle.Y[vecD]);
        lineMax = floor_441(triangle.Y[vecL]);
        //cout << "Down ";
    } else if (triangle.type == GoingUp) {
        lineMin = ceil_441(triangle.Y[vecL]);
        lineMax = floor_441(triangle.Y[vecD]);
        //cout << "Up ";
    }

    double leftSlope = scanline_slope(triangle.X[vecD], triangle.Y[vecD],
                                      triangle.X[vecL], triangle.Y[vecL]);
    double rightSlope = scanline_slope(triangle.X[vecD], triangle.Y[vecD],
                                       triangle.X[vecR], triangle.Y[vecR]);
    double leftB = scanline_b(triangle.Y[vecD], leftSlope, triangle.X[vecD]);
    double rightB = scanline_b(triangle.Y[vecD], rightSlope, triangle.X[vecD]);

    /*
    cout << "triangle: " << std::setprecision(16) << TRIANGLE_NUM << endl
         << "det: (" << triangle.X[vecD] << ", " << triangle.Y[vecD] << ")" << endl
         << "left: (" << triangle.X[vecL] << ", " << triangle.Y[vecL] << ")" << endl
         << "right: (" << triangle.X[vecR] << ", " << triangle.Y[vecR] << ")" << endl 
         << "r: " << (int) triangle.color[0] << ", g: " << (int) triangle.color[1] << ", b: " << (int) triangle.color[2] << endl
         << "lineMin = " << lineMin << ", lineMax = " << lineMax << endl
         << "l slope : " << leftSlope << ", r slope : " << std::setprecision(16) << rightSlope << endl
         << "l y-intercept : " << leftB << " r y-interept : " << rightB << endl << endl;
    
    */

    for (double line = lineMin; line <= lineMax; line++) {
        double leftEnd = ceil_441(scanline_end(triangle.X[vecD], line, leftSlope, leftB));
        double rightEnd = floor_441(scanline_end(triangle.X[vecD], line, rightSlope, rightB));

        for (double col = leftEnd; col <= rightEnd; col++)
            screen.setPixel(col, line, 0, triangle.color);
    }
}

void rasterizeArbitraryTriangle(Triangle triangle, Screen screen) {
    int vecD = triangle.vecD;
    int vecL = triangle.vecL;
    int vecR = triangle.vecR;

    Triangle upHalf = Triangle();
    Triangle downHalf = Triangle();

    double longSlope = scanline_slope(triangle.X[vecL], triangle.Y[vecL],
                                      triangle.X[vecR], triangle.Y[vecR]);
    double longB = scanline_b(triangle.Y[vecL], longSlope, triangle.X[vecL]);
    double longEnd = scanline_end(triangle.X[vecL], triangle.Y[vecD], longSlope, longB);

    /*
    cout << "ARB triangle: " << std::setprecision(16) << TRIANGLE_NUM << endl
         << "det: (" << triangle.X[vecD] << ", " << triangle.Y[vecD] << ")" << endl
         << "left: (" << triangle.X[vecL] << ", " << triangle.Y[vecL] << ")" << endl
         << "right: (" << triangle.X[vecR] << ", " << triangle.Y[vecR] << ")" << endl << endl; 
    */

    int newL = 0;
    int newR = 1;
    if (triangle.X[vecD] < longEnd) { newL = 1; newR = 0; }

    upHalf.X[0] = longEnd;
    upHalf.X[1] = triangle.X[vecD];
    upHalf.X[2] = triangle.X[vecR];
    upHalf.Y[0] = triangle.Y[vecD];
    upHalf.Y[1] = triangle.Y[vecD];
    upHalf.Y[2] = triangle.Y[vecR];
    upHalf.vecL = newL;
    upHalf.vecR = newR;
    upHalf.vecD = 2;
    upHalf.type = GoingUp;
    upHalf.color[0] = triangle.color[0];
    upHalf.color[1] = triangle.color[1];
    upHalf.color[2] = triangle.color[2];
    downHalf.X[0] = longEnd;
    downHalf.X[1] = triangle.X[vecD];
    downHalf.X[2] = triangle.X[vecL];
    downHalf.Y[0] = triangle.Y[vecD];
    downHalf.Y[1] = triangle.Y[vecD];
    downHalf.Y[2] = triangle.Y[vecL];
    downHalf.vecL = newL;
    downHalf.vecR = newR;
    downHalf.vecD = 2;
    downHalf.type = GoingDown;
    downHalf.color[0] = triangle.color[0];
    downHalf.color[1] = triangle.color[1];
    downHalf.color[2] = triangle.color[2];

    rasterizeTriangle(upHalf, screen);
    rasterizeTriangle(downHalf, screen);
}

void scanline(std::vector<Triangle> triangles, Screen screen) {
    for (auto & t : triangles) {
        setTriangleProps(&t);

        if (t.type == Arbitrary) { rasterizeArbitraryTriangle(t, screen); }
        else { rasterizeTriangle(t, screen); }

        TRIANGLE_NUM++;
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

    scanline(triangles, screen);

    WriteImage(image, "allTriangles");
}
