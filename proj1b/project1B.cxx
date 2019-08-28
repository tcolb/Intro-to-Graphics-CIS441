#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

const int WINDOW_HEIGHT = 1000;
const int WINDOW_WIDTH = 1000;

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

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      vtkImageData *image;
      unsigned char   *buffer;
      int width, height;

      // would some methods for accessing and setting pixels be helpful?
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
    std::vector<Triangle> rv(100);

    unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
    for (int i = 0 ; i < 100 ; i++)
    {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ;
       if (i == 5)
          rv[i].Y[(firstPt+2)%3] = -50;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
    }

    return rv;
}

double scanline_slope(double x1, double y1, double x2, double y2) {
    if (x1 == x2)
        return 0;
    return double((y2-y1)/(x2-x1));
}

double scanline_b(double y, double m, double x) {
    return double(y-m*x);
}


double scanline_line_min(double y) {
    return ceil_441(y);
}

double scanline_line_max(double y) {
    return floor_441(y);
}

double scanline_end(double x, double y, double m, double b) {
    if (m == 0)
        return x;
    return double((y-b)/m);
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

    //
    // begin scanline implementation
    //
    int triangleNum = 1; 
    for (auto & triangle : triangles) {

        //
        // determining vertex position relative to triangle
        //
        int topLeftVert, topRightVert, bottomVert = 0;
        double smallYVal = std::numeric_limits<double>::max();
        double smallXVal = std::numeric_limits<double>::max();
        double largeXVal = std::numeric_limits<double>::min();

        for (int i = 0; i < 3; i++)
            if (triangle.Y[i] < smallYVal) { smallYVal = triangle.Y[i]; bottomVert = i; }
        for (int i = 0; i < 3; i++) {
            if (i == bottomVert)
                continue;
            if (triangle.X[i] < smallXVal) { smallXVal = triangle.X[i]; topLeftVert = i; }
            if (triangle.X[i] > largeXVal) { largeXVal = triangle.X[i]; topRightVert = i; }
        }

        //
        // calclulating scanline values
        //
        double lineMin = scanline_line_min(triangle.Y[bottomVert]);
        double lineMax = scanline_line_max(triangle.Y[topLeftVert]);

        double leftSlope = scanline_slope(triangle.X[bottomVert], triangle.Y[bottomVert],
                                          triangle.X[topLeftVert], triangle.Y[topLeftVert]);
        double rightSlope = scanline_slope(triangle.X[bottomVert], triangle.Y[bottomVert],
                                           triangle.X[topRightVert], triangle.Y[topRightVert]);
        double leftB = scanline_b(triangle.Y[bottomVert], leftSlope, triangle.X[bottomVert]);
        double rightB = scanline_b(triangle.Y[bottomVert], rightSlope, triangle.X[bottomVert]);
       
        //
        // printing triangle data
        //
        cout << "triangle: " << std::setprecision(16) << triangleNum << endl
             << "bottom: (" << triangle.X[bottomVert] << ", " << triangle.Y[bottomVert] << ")" << endl
             << "left: (" << triangle.X[topLeftVert] << ", " << triangle.Y[topLeftVert] << ")" << endl
             << "right: (" << triangle.X[topRightVert] << ", " << triangle.Y[topRightVert] << ")" << endl 
             << "lineMin = " << lineMin << ", lineMax = " << lineMax << endl
             << "l slope : " << leftSlope << ", r slope : " << std::setprecision(16) << rightSlope << endl
             << "l y-intercept : " << leftB << " r y-interept : " << rightB << endl << endl;
    
        //
        // running scanline
        // 
        for (double line = lineMin; line <= lineMax; line++) {

            double leftEnd = ceil_441(scanline_end(triangle.X[bottomVert], line, leftSlope, leftB));
            double rightEnd = floor_441(scanline_end(triangle.X[bottomVert], line, rightSlope, rightB));
            
            for (double col = leftEnd; col <= rightEnd; col++) {
                screen.setPixel(col, line, 0, triangle.color);
                /*
                 screen.buffer = (unsigned char*) image->GetScalarPointer(col, line, 0);
                 screen.buffer[0] = triangle.color[0];
                 screen.buffer[1] = triangle.color[1];
                 screen.buffer[2] = triangle.color[2];
                 */
            }
        }

        triangleNum++;
    }

    WriteImage(image, "allTriangles");
}
