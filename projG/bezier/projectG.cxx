/**
 * CIS 541 Winter 2019
 * Project G
 * He He
 */

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
#include <vtkDoubleArray.h>
#include <cmath>

using std::cerr;
using std::endl;

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

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    void normalize_3vector(double* v)
    {
        double l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] = v[0]/l;
        v[1] = v[1]/l;
        v[2] = v[2]/l;
    } 

    // use your own project 1 code for those three functions
    Matrix          ViewTransform(void){
        Matrix view_tf;
        view_tf.A[0][0] = 1/tan(this->angle/2); view_tf.A[0][1] = 0;                    view_tf.A[0][2] = 0;                                               view_tf.A[0][3] = 0;
        view_tf.A[1][0] = 0;                    view_tf.A[1][1] = 1/tan(this->angle/2); view_tf.A[1][2] = 0;                                               view_tf.A[1][3] = 0;
        view_tf.A[2][0] = 0;                    view_tf.A[2][1] = 0;                    view_tf.A[2][2] = (this->far+this->near)/(this->far-this->near);   view_tf.A[2][3] = -1;
        view_tf.A[3][0] = 0;                    view_tf.A[3][1] = 0;                    view_tf.A[3][2] = (2*this->far*this->near)/(this->far-this->near); view_tf.A[3][3] = 0;
        return view_tf;
    };

    Matrix          CameraTransform(void){
        double O[3] = { this->position[0], this->position[1], this->position[2] };
        double T[3] = {-O[0], -O[1], -O[2]};
        double W[3] = {O[0]-this->focus[0], O[1]-this->focus[1], O[2]-this->focus[2]};
        double U[3] = {this->up[1]*W[2] - this->up[2]*W[1],
                       W[0]*this->up[2] - this->up[0]*W[2],
                       this->up[0]*W[1] - this->up[1]*W[0]};
        double V[3] = {W[1]*U[2] - W[2]*U[1],
                       U[0]*W[2] - W[0]*U[2],
                       W[0]*U[1] - W[1]*U[0]};

        normalize_3vector(U);
        normalize_3vector(V);
        normalize_3vector(W);

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
    };

    Matrix          DeviceTransform(int n, int m) {
        Matrix device_tf;
        device_tf.A[0][0] =  m/2; device_tf.A[0][1] =  0; device_tf.A[0][2] =  0; device_tf.A[0][3] = 0;
        device_tf.A[1][0] =  0; device_tf.A[1][1] =  n/2; device_tf.A[1][2] =  0; device_tf.A[1][3] = 0;
        device_tf.A[2][0] =  0; device_tf.A[2][1] =  0; device_tf.A[2][2] =  1; device_tf.A[2][3] = 0;
        device_tf.A[3][0] = m/2; device_tf.A[3][1] = n/2; device_tf.A[3][2] = 0; device_tf.A[3][3] = 1;
        return device_tf;
    };
};

// note that this is different from project 1 starter code
Camera
GetCamera(int frame, int nframes)
{
    double t = (double)frame/(double)nframes;
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 80 * sin(2*M_PI*t);
    c.position[1] = 30;
    c.position[2] = 80 * cos(2*M_PI*t);
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

class Screen
{
    public:
        vtkImageData    *image;
        unsigned char   *buffer;
        double          *zbuffer;
        int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

void colorPixel(Screen * screen, int y, int x, unsigned char color[3], double z){
    // implement this function, or use your own project 1E code.
    if (y < 0 || y >= screen->height || x < 0 || x >= screen->width)
        return;

    
    if (screen->zbuffer[y*x] > z)
      return;

    screen->zbuffer[y*x] = z;
    screen->buffer = (unsigned char*) screen->image->GetScalarPointer(x, y , 0);
    screen->buffer[0] = color[0];
    screen->buffer[1] = color[1];
    screen->buffer[2] = color[2];

}

class Line{
public:
    double X[2];
    double Y[2];
    double Z[2];
    unsigned char   color[3];

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }
    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Draw(Screen * screen){
        int x0, y0, x1, y1;
        if(Y[0] < Y[1]){
            x0 = X[0];
            y0 = Y[0];
            x1 = X[1];
            y1 = Y[1];
        }else{
            x0 = X[1];
            y0 = Y[1];
            x1 = X[0];
            y1 = Y[0];
        }

        int e;
        int dx = x1 - x0;
        int dy = y1 - y0;
        int x = x0;
        int y = y0;
        double dz = Z[1] - Z[0];

        int stepx = 1;
        int stepy = 1;
        if (dx<0){
            dx = - dx;
            stepx = -1;
        }
        if (dy<0){
            dy = - dy;
            stepy = -1;
        }

        int i;
        // first octant
        if(dy<dx){
            e = 2 * dy - dx;
            for(i=0; i<dx; i++){
                colorPixel(screen, y, x, color, i*dz/dx+Z[0]);
                while(e > 0){
                    y += stepy;
                    e = e - 2 * dx;
                }
                x += stepx;
                e = e + 2 * dy;
            }
        // second octant
        }else{
            e = 2 * dx - dy;
            for(i=0; i<dy; i++){
                colorPixel(screen, y, x, color, i*dz/dy+Z[0]);
                while(e > 0){
                    x += stepx;
                    e = e - 2 * dy;
                }
                y += stepy;
                e = e + 2 * dx;
            }
        }
    }

    void Print(){
        std::cout << "v1: " << X[0] << ", " << Y[0] << ", " << Z[0] << '\n';
        std::cout << "v2: " << X[1] << ", " << Y[1] << ", " << Z[1] << '\n';
    }
};

void casteljau_helper(double p[4], double q[4], double r[4]) {
    q[0] = p[0];
    q[1] = 0.5*(p[0] + p[1]);
    q[2] = 0.5*q[1] + 0.25*(p[1] + p[2]);
    r[3] = p[3];
    r[2] = 0.5*(p[2] + p[3]);
    r[1] = 0.5*r[2] + 0.25*(p[1] + p[2]);
    q[3] = 0.5*(q[2] + r[1]);
    r[0] = q[3];
}

void Bezier_Divide(double pX[4], double pY[4], double pZ[4], double qX[4], double qY[4], double qZ[4], double rX[4], double rY[4], double rZ[4]){
    // your code goes here
    casteljau_helper(pX, qX, rX);
    casteljau_helper(pY, qY, rY);
    casteljau_helper(pZ, qZ, rZ);
}

class BezierCurve{
public:
    double X[4];
    double Y[4];
    double Z[4];
    unsigned char color[3];

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Print(){
        printf("p0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("p1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("p2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);
        printf("p3: %lf, %lf, %lf\n", X[3], Y[3], Z[3]);
    }

    void RotateY(BezierCurve * dst, double angleDegree){
        int i;
        for(i=0; i<4; i++){
            dst->X[i] = X[i] * cos(2*M_PI*angleDegree/360)
                        + Z[i] * sin(2*M_PI*angleDegree/360);
            dst->Y[i] = Y[i];
            dst->Z[i] = Z[i] * cos(2*M_PI*angleDegree/360)
                        - X[i] * sin(2*M_PI*angleDegree/360);
        }
    }

    void Draw(Screen * screen){
        // your code goes here
        int Xmin, Xmax, Ymin, Ymax;
        Xmin = screen->width;
        Ymin = screen->height;
        Xmax = Ymax = 0;

        // determine mins / maxs
        for (int i = 0; i < 4; i++) {
            if (X[i] > Xmax)
                Xmax = X[i];
            if (X[i] < Xmin)
                Xmin = X[i];
            if (Y[i] > Ymax)
                Ymax = X[i];
            if (Y[i] < Ymin)
                Ymin = X[i];
        }

        if (Ymax - Ymin < 10 && Xmax - Xmin < 10) {
            // do drawing
            for (int i = 0; i < 3; i++) {
                Line l;
                l.X[0] = this->X[i];
                l.X[1] = this->X[i+1];
                l.Y[0] = this->Y[i];
                l.Y[1] = this->Y[i+1];
                l.Z[0] = this->Z[i];
                l.Z[1] = this->Z[i+1];
                l.color[0] = this->color[0];
                l.color[1] = this->color[1];
                l.color[2] = this->color[2];
                //cout << "line : " << "(" << l.X[0] << ", " << l.Y[0] << ", " << l.Z[0] << "), ("<< l.X[1] << ", " << l.Y[1] << ", " << l.Z[1] << ")" << endl;

                l.Draw(screen);
            }
        } else {
            // do splitting
            double qX[4], qY[4], qZ[4], rX[4], rY[4], rZ[4];
            Bezier_Divide(this->X, this->Y, this->Z, qX, qY, qZ, rX, rY, rZ);
            BezierCurve first;
            BezierCurve second;

            for (int i = 0; i < 4; i++) {
                first.X[i] = qX[i];
                first.Y[i] = qY[i];
                first.Z[i] = qZ[i];
                second.X[i] = rX[i];
                second.Y[i] = rY[i];
                second.Z[i] = rZ[i];
            }

            for (int i = 0; i < 3; i++) {
                first.color[i] = this->color[i];
                second.color[i] = this->color[i];
            }

            first.Draw(screen);
            second.Draw(screen);
        }

        return;
    }
};

class BezierSurface{
public:
    double X[16];
    double Y[16];
    double Z[16];
    unsigned char color[3];

    void SetPoints(double x[16], double y[16], double z[16]){
        int i;
        for(i=0; i<16; i++){
            X[i] = x[i];
            Y[i] = y[i];
            Z[i] = z[i];
        }
    }

    void SetPointByIndex(int index, double x, double y, double z){
        X[index] = x;
        Y[index] = y;
        Z[index] = z;
    }

    void SetColor(unsigned char r, unsigned char g, unsigned char b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void Print(){
        printf("v0: %lf, %lf, %lf\n", X[0], Y[0], Z[0]);
        printf("v1: %lf, %lf, %lf\n", X[1], Y[1], Z[1]);
        printf("v2: %lf, %lf, %lf\n", X[2], Y[2], Z[2]);

    }

    void Draw(Screen * screen, int divisions){
        // your code goes here
        // do divisons
        if (divisions > 0) {
            BezierSurface tlSurface, trSurface, blSurface, brSurface, lIntermediate, rIntermediate;
            // calc divisons
            double qX[4], qY[4], qZ[4], rX[4], rY[4], rZ[4];
            for (int r=0; r<4; r++) {
                int idx = 4*r;
                Bezier_Divide(this->X + idx, this->Y + idx, this->Z + idx, qX, qY, qZ, rX, rY, rZ);
                for (int i=0; i<4; i++) {
                    lIntermediate.X[idx + i] = qX[i];
                    lIntermediate.Y[idx + i] = qY[i];
                    lIntermediate.Z[idx + i] = qZ[i];
                    rIntermediate.X[idx + i] = rX[i];
                    rIntermediate.Y[idx + i] = rY[i];
                    rIntermediate.Z[idx + i] = rZ[i];
                }
            }

            for (int c=0; c<4; c++) {
                int idx;
                double tX[4], tY[4], tZ[4];
                for (int i=0; i<4; i++) {
                    idx = 4*i;
                    tX[i] = lIntermediate.X[idx + c];
                    tY[i] = lIntermediate.Y[idx + c];
                    tZ[i] = lIntermediate.Z[idx + c];
                }

                Bezier_Divide(tX, tY, tZ, qX, qY, qZ, rX, rY, rZ);
                for (int i=0; i<4; i++) {
                    idx = 4*i;
                    tlSurface.X[idx + c] = qX[i];
                    tlSurface.Y[idx + c] = qY[i];
                    tlSurface.Z[idx + c] = qZ[i];
                    blSurface.X[idx + c] = rX[i];
                    blSurface.Y[idx + c] = rY[i];
                    blSurface.Z[idx + c] = rZ[i];
                }

                for (int i=0; i<4; i++) {
                    idx = 4*i;
                    tX[i] = rIntermediate.X[idx + c];
                    tY[i] = rIntermediate.Y[idx + c];
                    tZ[i] = rIntermediate.Z[idx + c];
                }

                Bezier_Divide(tX, tY, tZ, qX, qY, qZ, rX, rY, rZ);
                for (int i=0; i<4; i++) {
                    idx = 4*i;
                    trSurface.X[idx + c] = qX[i];
                    trSurface.Y[idx + c] = qY[i];
                    trSurface.Z[idx + c] = qZ[i];
                    brSurface.X[idx + c] = rX[i];
                    brSurface.Y[idx + c] = rY[i];
                    brSurface.Z[idx + c] = rZ[i];
                }
            }

            tlSurface.SetColor(this->color[0],  this->color[1], this->color[2]);
            trSurface.SetColor(this->color[0],  this->color[1], this->color[2]);
            blSurface.SetColor(this->color[0],  this->color[1], this->color[2]);
            brSurface.SetColor(this->color[0],  this->color[1], this->color[2]);

            // recurse divisions
            tlSurface.Draw(screen, divisions-1);
            trSurface.Draw(screen, divisions-1);
            blSurface.Draw(screen, divisions-1);
            brSurface.Draw(screen, divisions-1);
        } else {
            // draw
            for (int i=0; i<4; i++) {
                // draw rows
                int idx = i*4;
                for (int r=0; r<3; r++) {
                    Line l;
                    l.X[0] = this->X[idx + r];
                    l.X[1] = this->X[idx + r+1];
                    l.Y[0] = this->Y[idx + r];
                    l.Y[1] = this->Y[idx + r+1];
                    l.Z[0] = this->Z[idx + r];
                    l.Z[1] = this->Z[idx + r+1];
                    l.color[0] = this->color[0];
                    l.color[1] = this->color[1];
                    l.color[2] = this->color[2];
                    l.Draw(screen);
                }
                
                // draw cols
                for (int c=0; c<3; c++) {
                    idx = c * 4;
                    Line l;
                    l.X[0] = this->X[idx + i];
                    l.X[1] = this->X[idx+4 + i];
                    l.Y[0] = this->Y[idx + i];
                    l.Y[1] = this->Y[idx+4 + i];
                    l.Z[0] = this->Z[idx + i];
                    l.Z[1] = this->Z[idx+4 + i];
                    l.color[0] = this->color[0];
                    l.color[1] = this->color[1];
                    l.color[2] = this->color[2];
                    l.Draw(screen);
                }
            }
        }

        return;
    }

};

// returns an array of BezierCurves of size 6
BezierCurve * getLeaves(){
    BezierCurve *c = (BezierCurve*)malloc(sizeof(BezierCurve)*6);
    int i;
    for(i=0; i<6; i++)
        c[i].SetColor(30, 215, 97);

    c[0].SetPointByIndex(0, 0, -10, 0);
    c[0].SetPointByIndex(1, 14, -7.2, 0);
    c[0].SetPointByIndex(2, 9.8, -3, 2.8);
    c[0].SetPointByIndex(3, 14, 4, 14);

    c[1].SetPointByIndex(0, 0, -10, 0);
    c[1].SetPointByIndex(1, 0, -7.2, 14);
    c[1].SetPointByIndex(2, 2.8, -3, 9.8);
    c[1].SetPointByIndex(3, 14, 4, 14);

    c[0].RotateY(&(c[2]), 120);
    c[1].RotateY(&(c[3]), 120);
    c[0].RotateY(&(c[4]), 240);
    c[1].RotateY(&(c[5]), 240);

    return &(c[0]);
}

// returns an array of BezierSurfaces of size 2
BezierSurface * getSurfaces(){
    BezierSurface *s = (BezierSurface*)malloc(sizeof(BezierSurface)*2);

    s[0].SetPointByIndex(0, 0.0, 0.0, 0.0); // first row
	s[0].SetPointByIndex(1, 0.0, 3, 5);
	s[0].SetPointByIndex(2, 0.0, 7.5, 10);
	s[0].SetPointByIndex(3, 0.0, 1.5, 15.0);
	s[0].SetPointByIndex(4, 5, 12, 0.0); // second row
	s[0].SetPointByIndex(5, 5, -1.5, 5);
	s[0].SetPointByIndex(6, 5, 0.0, 10);
	s[0].SetPointByIndex(7, 5, 1.5, 15.0);
	s[0].SetPointByIndex(8, 10, 4.5, 0.0); // third row
	s[0].SetPointByIndex(9, 10, 12, 5);
	s[0].SetPointByIndex(10, 10, 13.5, 10);
	s[0].SetPointByIndex(11, 10, 7.5, 15.0);
	s[0].SetPointByIndex(12, 15.0, 6, 0.0); // fourth row
	s[0].SetPointByIndex(13, 15.0, 3, 5);
	s[0].SetPointByIndex(14, 15.0, 7.5, 10);
	s[0].SetPointByIndex(15, 15.0, 15.0, 15.0);
    s[0].SetColor(51, 133, 229);

    s[1].SetPointByIndex(0, 0.0, -3, 0.0); // first row
	s[1].SetPointByIndex(1, 0.0, -3, 5);
	s[1].SetPointByIndex(2, 0.0, -3, 10);
	s[1].SetPointByIndex(3, 0.0, -3, 15);
	s[1].SetPointByIndex(4, 5, -3, 0.0); // second row
	s[1].SetPointByIndex(5, 5, -3, 5);
	s[1].SetPointByIndex(6, 5, -3, 10);
	s[1].SetPointByIndex(7, 5, -3, 15);
	s[1].SetPointByIndex(8, 10, -3, 0.0); // third row
	s[1].SetPointByIndex(9, 10, -3, 5);
	s[1].SetPointByIndex(10, 10, -3, 10);
	s[1].SetPointByIndex(11, 10, -3, 15);
	s[1].SetPointByIndex(12, 15, -3, 0.0); // fourth row
	s[1].SetPointByIndex(13, 15, -3, 5);
	s[1].SetPointByIndex(14, 15, -3, 10);
	s[1].SetPointByIndex(15, 15, -3, 15);
    s[1].SetColor(31, 179, 83);

    return &(s[0]);
}

BezierSurface *getTeapotSurfaces() {
    BezierSurface *t = (BezierSurface*)malloc(sizeof(BezierSurface)*32);
    for (int i=0; i<32; i++) {
        double buf = 0;
        cin >> buf;
        cin >> buf;
        for (int j=0; j<16; j++) {
            cin >> buf;
            t[i].X[j] = buf*5;
            cin >> buf;
            t[i].Y[j] = buf*5-3.5;
            cin >> buf;
            t[i].Z[j] = buf*5;
        }
        t[i].SetColor(255, 255, 255);
    }
    return &(t[0]);
}

int main(int argc, char** argv)
{
    int choice = 0;
    if (argc == 2) {
        choice = atoi(argv[1]);
    }

    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer =
        (unsigned char *) image->GetScalarPointer(0,0,0);

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    screen.image = image;
    screen.zbuffer = (double*)malloc(sizeof(double) * screen.width * screen.height);

    // uncooment the following two lines to get the curves and surfaces
    BezierCurve *s = getLeaves();
    BezierSurface *b;
    if (choice == 2) {
        b = getTeapotSurfaces();
    } else {
        b = getSurfaces();
    }
    
    // Camera rotates around y-axis
    // It comes back to the original position after 100 iterations
    int im;
    for (im=0; im<100; im++) {
        cout << "setting image " << im << endl;
        int npixels = screen.width * screen.height;
        screen.buffer = (unsigned char*) screen.image->GetScalarPointer(0,0,0);
        for (int j = 0 ; j < npixels*3 ; j++)
            screen.buffer[j] = 0;
        for(int j = 0; j < npixels; j++)
            screen.zbuffer[j] = -1;

        Camera camera = GetCamera(im*10, 1000);
        Matrix ct = camera.CameraTransform();
        Matrix vt = camera.ViewTransform();
        Matrix intermediate = Matrix::ComposeMatrices(ct, vt);
        Matrix dt = camera.DeviceTransform(screen.width, screen.height);
        Matrix total = Matrix::ComposeMatrices(intermediate, dt);

        switch(choice) {
            case 0: {

        for (int i=0; i<6; i++) {
            //cout << "transforming curve " << i << endl;
            BezierCurve tfCurve;
            for (int j=0; j<4; j++) {
                double pt[4];
                double o_pt[4];
    
                pt[0] = s[i].X[j];
                pt[1] = s[i].Y[j];
                pt[2] = s[i].Z[j];
                pt[3] = 1;

                total.TransformPoint(pt, o_pt);


                tfCurve.X[j] = o_pt[0]/o_pt[3];
                tfCurve.Y[j] = o_pt[1]/o_pt[3];
                tfCurve.Z[j] = o_pt[2]/o_pt[3];
                tfCurve.color[j] = s[i].color[j];
                //std::cout << tfCurve.X[j] << ", " << tfCurve.Y[j] << ", " << tfCurve.Z[j] << std::endl;
            }
            tfCurve.Draw(&screen);
        }
        break;
            }

            case 1: {
        for (int i = 0; i < 2; i++) {
            BezierSurface tfSurface;
            for (int j=0; j<16; j++) {
                double pt[4];
                double o_pt[4];
                pt[0] = b[i].X[j];
                pt[1] = b[i].Y[j];
                pt[2] = b[i].Z[j];
                pt[3] = 1;

                total.TransformPoint(pt, o_pt);

                tfSurface.X[j] = o_pt[0]/o_pt[3];
                tfSurface.Y[j] = o_pt[1]/o_pt[3];
                tfSurface.Z[j] = o_pt[2]/o_pt[3];
                //std::cout << "in pt: " << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
                //std::cout << "out pt: " << tfSurface.X[j] << ", " << tfSurface.Y[j] << ", " << tfSurface.Z[j] << std::endl;
            }
            tfSurface.SetColor(b[i].color[0], b[i].color[1], b[i].color[2]);
            tfSurface.Draw(&screen, 3);
        }
        break;
            }
            case 2: {
       // cout << screen.height << " " << screen.width << endl;
       for (int i=0; i<32; i++) {
            BezierSurface tfSurface;
            for (int j=0; j<16; j++) {
                double pt[4];
                double o_pt[4];
                pt[0] = b[i].X[j];
                pt[1] = b[i].Y[j];
                pt[2] = b[i].Z[j];
                pt[3] = 1;

                total.TransformPoint(pt, o_pt);

                tfSurface.X[j] = o_pt[0]/o_pt[3];
                tfSurface.Y[j] = o_pt[1]/o_pt[3];
                tfSurface.Z[j] = o_pt[2]/o_pt[3];
            }
            tfSurface.SetColor(b[i].color[0], b[i].color[1], b[i].color[2]);
            tfSurface.Draw(&screen, 3);
       }
       break;
            }
        }

        char name[20];
        // make sure you have a directory named images
        // so that writing images won't cause you errors
        sprintf(name, "./images/frame%02d", im);
        WriteImage(image, name);
    }
}
