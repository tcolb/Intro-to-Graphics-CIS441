#include <iostream>
#include <array>
#include <cmath>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;


enum TriangleType {GoingUp, GoingDown, Arbitrary};
typedef std::array<double, 3> tvec;
const int WINDOW_HEIGHT = 1000;
const int WINDOW_WIDTH = 1000;
int TRIANGLE_NUM = 0;


//
// vtk helpers
//

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


//
// triangle and screen
//

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
          this->buffer = (unsigned char*) image->GetScalarPointer(x, y , 0);
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


//
// matrix and helpers
//

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};


void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}


Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}


void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}


//
// camera and helpers
//

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void)
    {
        Matrix view_tf;
        view_tf.A[0][0] = 1/tan(this->angle/2); view_tf.A[0][1] = 0;                view_tf.A[0][2] = 0; view_tf.A[0][3] = 0;
        view_tf.A[1][0] = 0;                view_tf.A[1][1] = 1/tan(this->angle/2); view_tf.A[1][2] = 0; view_tf.A[1][3] = 0;
        view_tf.A[2][0] = 0;                view_tf.A[2][1] = 0;                view_tf.A[2][2] = (this->far+this->near)/(this->far-this->near); view_tf.A[2][3] = -1;
        view_tf.A[3][0] = 0;                view_tf.A[3][1] = 0;                view_tf.A[3][2] = (2*this->far*this->near)/(this->far-this->near); view_tf.A[3][3] = 0;
        return view_tf;
    }

    Matrix          CameraTransform(void)
    {
        double O[3] = { this->position[0], this->position[1], this->position[2] };
        double T[3] = {-O[0], -O[1], -O[2]};
        double W[3] = {O[0]-this->focus[0], O[1]-this->focus[1], O[2]-this->focus[2]};
        double U[3] = {this->up[1]*W[2] - this->up[2]*W[1],
                       W[0]*this->up[2] - this->up[0]*W[2],
                       this->up[0]*W[1] - this->up[1]*W[0]};
        double V[3] = {W[1]*U[2] - W[2]*U[1],
                       U[0]*W[2] - W[0]*U[2],
                       W[0]*U[1] - W[1]*U[0]};

        normalize_vector(U);
        normalize_vector(V);
        normalize_vector(W);

        #ifdef DEBUG_TF
        cout << "Camera Frame: U = " << U[0] << ", " << U[1] << ", " << U[2] << endl;
        cout << "Camera Frame: V = " << V[0] << ", " << V[1] << ", " << V[2] << endl;
        cout << "Camera Frame: W = " << W[0] << ", " << W[1] << ", " << W[2] << endl;
        cout << "Camera Frame: O = " << O[0] << ", " << O[1] << ", " << O[2] << endl;
        #endif

        double UT[3] = {U[0]*T[0] + U[1]*T[1] + U[2]*T[2],
                        V[0]*T[0] + V[1]*T[1] + V[2]*T[2],
                        W[0]*T[0] + W[1]*T[1] + W[2]*T[2]};

        // camera transform
        Matrix camera_tf;
        camera_tf.A[0][0] =  U[0]; camera_tf.A[0][1] =  V[0]; camera_tf.A[0][2] =  W[0]; camera_tf.A[0][3] = 0;
        camera_tf.A[1][0] =  U[1]; camera_tf.A[1][1] =  V[1]; camera_tf.A[1][2] =  W[1]; camera_tf.A[1][3] = 0;
        camera_tf.A[2][0] =  U[2]; camera_tf.A[2][1] =  V[2]; camera_tf.A[2][2] =  W[2]; camera_tf.A[2][3] = 0;
        camera_tf.A[3][0] = UT[0]; camera_tf.A[3][1] = UT[1]; camera_tf.A[3][2] = UT[2]; camera_tf.A[3][3] = 1;
        return camera_tf;
    }

    Matrix          DeviceTransform(void)
    {
        Matrix device_tf;
        zeroMatrix(&device_tf);
        device_tf.A[0][0] = WINDOW_WIDTH/2;
        device_tf.A[1][1] = WINDOW_HEIGHT/2;
        device_tf.A[2][2] = 1;
        device_tf.A[3][0] = WINDOW_WIDTH/2;
        device_tf.A[3][1] = WINDOW_HEIGHT/2;
        device_tf.A[3][3] = 1;
        return device_tf;
    }

    Matrix          CompositeTransform(void)
    {
        return Matrix::ComposeMatrices(Matrix::ComposeMatrices(this->CameraTransform(), this->ViewTransform()), this->DeviceTransform());
    }

  private:
    void normalize_vector(double *v)
    {
        double l = sqrt(pow(v[0], 2.0) + pow(v[1], 2.0) + pow(v[2], 2.0));
        v[0] = v[0]/l;
        v[1] = v[1]/l;
        v[2] = v[2]/l;
    }

    void zeroMatrix(Matrix* m)
    {
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                m->A[i][j] = 0;
    }
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}


Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}


//
// generate triangles
//

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
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
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
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
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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

std::vector<Triangle>
TransformTriangles(std::vector<Triangle> tris, Camera c)
{
    std::vector<Triangle> outTris;
    Matrix M = c.CompositeTransform();
    for (auto & t : tris)
    {
	#ifdef DEBUG_T
	cout << "TRIANGLE: " << std::setprecision(16) << TRIANGLE_NUM << endl
	     << "BEFORE TRANSFORM:" << endl
	     << "det: (" << t.X[t.vecD] << ", " << t.Y[t.vecD] << ", " << t.Z[t.vecD] << ") "
	     << "left: (" << t.X[t.vecL] << ", " << t.Y[t.vecL] << ", " << t.Z[t.vecL] << ") "
	     << "right: (" << t.X[t.vecR] << ", " << t.Y[t.vecR] << ", " << t.Z[t.vecR] << ") " << endl << endl;
	#endif

        Triangle newT;
        newT.colors = t.colors;
        newT.type = t.type;
        newT.vecD = t.vecD;
        newT.vecL = t.vecL;
        newT.vecR = t.vecR;

	for (int i=0; i<3; i++)
        {
	    double I_vert[4] = { t.X[i], t.Y[i], t.Z[i], 1.0 };
	    double M_vert[4] = { 0 };
	    M.TransformPoint(I_vert, M_vert);
	    newT.X[i] = M_vert[0]/M_vert[3];
	    newT.Y[i] = M_vert[1]/M_vert[3];
	    newT.Z[i] = M_vert[2]/M_vert[3];
        }
        
	#ifdef DEBUG_T
	cout << "AFTER TRANSFORM:" << endl
	     << "denewT: (" << newT.X[newT.vecD] << ", " << newT.Y[newT.vecD] << ", " << newT.Z[newT.vecD] << ") "
	     << "lefnewT: (" << newT.X[newT.vecL] << ", " << newT.Y[newT.vecL] << ", " << newT.Z[newT.vecL] << ") "
	     << "righnewT: (" << newT.X[newT.vecR] << ", " << newT.Y[newT.vecR] << ", " << newT.Z[newT.vecR] << ") " << endl << endl;
	#endif

        outTris.push_back(newT);
    }
    return outTris;	
}

//
// general helpers
//

double calcT(double x, double a, double b)
{
    return double((x-a)/(b-a));   
}

double lerp(double t, double fa, double fb)
{
    return double(fa + t*(fb-fa));
}

void color_lerp(double t, tvec aColor, tvec bColor, tvec &out )
{
    out[0] = lerp(t, aColor[0], bColor[0]);
    out[1] = lerp(t, aColor[1], bColor[1]);
    out[2] = lerp(t, aColor[2], bColor[2]);
}

void arr_cpy(double* a1, double* a2, int n1, int n2)
{
    for (int i=0; i<n1; i++) {
        if (i > n2) return;
        a2[i] = a1[i];
    }
}

void setTriangleProps(Triangle* t)
{
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


//
// rasterization and scanline
//

void rasterizeTriangle(Triangle t, Screen* screen)
{
    // fix going up / down min and max
    double rowMin, rowMax;
    
    if (t.type == GoingDown) {
        rowMin = ceil_441(t.Y[t.vecD]);
        rowMax = floor_441(t.Y[t.vecL]);
        //cout << "Down ";
    } else if (t.type == GoingUp) {
        rowMin = ceil_441(t.Y[t.vecL]);
        rowMax = floor_441(t.Y[t.vecD]);
        //cout << "Up ";
    }
  
    // scanline row
    for (double row = rowMin; row <= rowMax; row++) {

        double lerpT = calcT(row, t.Y[t.vecD], t.Y[t.vecL]);   
        // not ceil/floored as using for color lerp as well
        double lEnd = lerp(lerpT, t.X[t.vecD], t.X[t.vecL]);
        double rEnd = lerp(lerpT, t.X[t.vecD], t.X[t.vecR]);
        // color lerping
        tvec lColor, rColor;
        color_lerp(lerpT, t.colors[t.vecD], t.colors[t.vecL], lColor);
        color_lerp(lerpT, t.colors[t.vecD], t.colors[t.vecR], rColor);
        // z-buffer lerping
        double lZ = lerp(lerpT, t.Z[t.vecD], t.Z[t.vecL]);
        double rZ = lerp(lerpT, t.Z[t.vecD], t.Z[t.vecR]);

        // scanline col
        for (double col = ceil_441(lEnd); col <= floor_441(rEnd); col++) {
            double xT = calcT(col, lEnd, rEnd);
            double xZ = lerp(xT, lZ, rZ);
            tvec xColor; 
            color_lerp(xT, lColor, rColor, xColor);
            screen->setPixel(col, row, xZ, xColor);
        }
    }
}

void rasterizeArbitraryTriangle(Triangle t, Screen* screen)
{
    Triangle upHalf = Triangle();
    Triangle downHalf = Triangle();

    double longT = calcT(t.Y[t.vecD], t.Y[t.vecL], t.Y[t.vecR]);
    double longEnd = lerp(longT, t.X[t.vecL], t.X[t.vecR]);

    #ifdef DEBUG_ARB
    cout << "ARB tangle: " << std::setprecision(16) << TRIANGLE_NUM << endl
         << "det: (" << t.X[t.vecD] << ", " << t.Y[t.vecD] << ", " << t.Z[t.vecD] << ")" << endl
         << "left: (" << t.X[t.vecL] << ", " << t.Y[t.vecL] << ", " << t.Z[t.vecL] << ")" << endl
         << "right: (" << t.X[t.vecR] << ", " << t.Y[t.vecR] << ", " << t.Z[t.vecR] << ")" << endl << endl; 
    #endif

    // fixing right / left verts
    int newL = 0;
    int newR = 1;
    if (t.X[t.vecD] < longEnd) { newL = 1; newR = 0; }

    // lerp endpoint colors for new ts
    tvec longColor;
    color_lerp(longT, t.colors[t.vecL], t.colors[t.vecR], longColor);

    // lerp endpoint z values for new ts
    double longZ = lerp(longT, t.Z[t.vecL], t.Z[t.vecR]);

    // setting going up
    upHalf.X[0] = longEnd;
    upHalf.X[1] = t.X[t.vecD];
    upHalf.X[2] = t.X[t.vecR];
    upHalf.Y[0] = t.Y[t.vecD];
    upHalf.Y[1] = t.Y[t.vecD];
    upHalf.Y[2] = t.Y[t.vecR];
    upHalf.Z[0] = longZ;
    upHalf.Z[1] = t.Z[t.vecD];
    upHalf.Z[2] = t.Z[t.vecR];
    upHalf.vecL = newL;
    upHalf.vecR = newR;
    upHalf.vecD = 2;
    upHalf.type = GoingUp;
    upHalf.colors[0] = longColor;
    upHalf.colors[1] = t.colors[t.vecD];
    upHalf.colors[2] = t.colors[t.vecR];

    // setting going down
    downHalf.X[0] = longEnd;
    downHalf.X[1] = t.X[t.vecD];
    downHalf.X[2] = t.X[t.vecL];
    downHalf.Y[0] = t.Y[t.vecD];
    downHalf.Y[1] = t.Y[t.vecD];
    downHalf.Y[2] = t.Y[t.vecL];
    downHalf.Z[0] = longZ;
    downHalf.Z[1] = t.Z[t.vecD];
    downHalf.Z[2] = t.Z[t.vecL];
    downHalf.vecL = newL;
    downHalf.vecR = newR;
    downHalf.vecD = 2;
    downHalf.type = GoingDown;
    downHalf.colors[0] = longColor;
    downHalf.colors[1] = t.colors[t.vecD];
    downHalf.colors[2] = t.colors[t.vecL];

    rasterizeTriangle(upHalf, screen);
    rasterizeTriangle(downHalf, screen);
}

void scanline(std::vector<Triangle> triangles, Screen* screen)
{
    for (auto & t : triangles) {
        setTriangleProps(&t);

        if (t.type == Arbitrary) { rasterizeArbitraryTriangle(t, screen); }
        else { rasterizeTriangle(t, screen); }

        TRIANGLE_NUM++;
    }
}


//
// mainline
//

int main()
{
    vtkImageData *image = NewImage(WINDOW_WIDTH, WINDOW_HEIGHT);
    unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = WINDOW_WIDTH*WINDOW_HEIGHT;
    for (int i = 0 ; i < npixels*3 ; i++) // make all pixels black
       buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();
 
    Screen screen;
    screen.image = image;
    screen.buffer = buffer;
    screen.width = WINDOW_WIDTH;
    screen.height = WINDOW_HEIGHT;


    for (int i=0; i<4; i++)
    {
        Camera c = GetCamera(i*250,1000);
        std::vector<Triangle> deviceTriangles = TransformTriangles(triangles, c); // apply composite transform vector on triangles
    	scanline(deviceTriangles, &screen);

        char filename[32];
        sprintf(filename, "frame%03d", i*250);
    	WriteImage(image, filename);
    }
        

    #ifdef DEBUG_TF
    Camera camera = GetCamera(0,1000);
    std::vector<Triangle> deviceTriangles = TransformTriangles(triangles, camera); // apply composite transform vector on triangles
    cout << "angle: " << camera.angle << ", far: " << camera.far << ", near: " << camera.near << endl;
    cout << "View transform matrix" << endl;
    view_tf.Print(cout);
    cout << endl;
    cout << "Device transform matrix" << endl;
    device_tf.Print(cout);
    cout << endl; 
    cout << "Camera transform matrix" << endl;
    camera_tf.Print(cout);
    cout << endl;
    cout << "Composite transform matrix" << endl;
    M_tf.Print(cout);
    cout << endl;
    

    double I_vert[4] = { 0, 36.4645, 36.4645, 1 };
    double M_vert[4] = { 0 };
    
    M_tf.TransformPoint(I_vert, M_vert);
    M_vert[0] = M_vert[0]/M_vert[3];
    M_vert[1] = M_vert[1]/M_vert[3];
    M_vert[2] = M_vert[2]/M_vert[3];
    cout << "Test Transform = " << M_vert[0] << ", " << M_vert[1] << ", " << M_vert[2] << endl;
    #endif

}

