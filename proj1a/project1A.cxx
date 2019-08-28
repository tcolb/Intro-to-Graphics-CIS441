#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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


int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350);
   int imgDims[3];
   image->GetDimensions(imgDims);
   int modRange[] = {0, 128, 255};
   int stripHeight = imgDims[1] / 27;

   //std::cerr << "Dims: " << imgDims[0] << ", " << imgDims[1] << " StripHeight: " << stripHeight << endl;
   
   int stripNum = 0;
   for (int y = 0; y < imgDims[1]; y = y + stripHeight) {
   	unsigned char *buffer;
	int b = modRange[stripNum%3];
	int g = modRange[(stripNum/3)%3];
	int r = modRange[stripNum/9];
   	for (int i = 0; i < stripHeight; i++) {
		//std::cerr << "strip " << stripNum << " line " << y+i << ", " << r << " " << g << " " << b << endl; 
		for (int x = 0; x < imgDims[0]; x++) {
			buffer = (unsigned char *) image->GetScalarPointer(x,y+i,0);
			buffer[0] = r; buffer[1] = g; buffer[2] = b;
		}
	}
	stripNum++;
   }

   WriteImage(image, "proj1A");
}
