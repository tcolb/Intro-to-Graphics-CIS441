#include <iostream>
#include <string>
#include <cstdio>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <cuda_runtime.h>
#include <iomanip>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cassert>
#include <cmath>

__global__
void rgba_to_greyscale(uchar4* d_rgba, unsigned char* d_grey, int N)
{
    // Don't forget to check if the index is out of bounds
    // A simple `return` will break out for us

    // Suggest you use a static_cast when converting back to your
    // grey image's index

    // L = 0.21 ∗ Red + 0.72 ∗ Green + 0.07 ∗ Blue
    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < N) {
       d_grey[idx] = 0.21 * d_rgba[idx].x + 0.72 * d_rgba[idx].y + 0.07 * d_rgba[idx].z; 
    }
    	
    return;
}

__global__ void
sobel_filter(unsigned char* d_grey, unsigned char* d_sobel, int rows, int cols)
{
    extern __shared__ uchar s_input2[];

    int sobelWidth = 3;
    int dx[3][3] = {-1, 0, 1,
                    -2, 0, 2,
                    -1, 0, 1};
    int dy[3][3] = { 1, 2, 1,
                     0, 0, 0,
                    -1,-2,-1};
    int sum_x = 0;
    int sum_y = 0;
    const int r = (sobelWidth - 1)/2;
    
    int block_size = blockDim.x*blockDim.y;
    int block_y = ((block_size)*gridDim.x)*blockIdx.y;
    int block_x = blockIdx.x*blockDim.x;
    int thread_y = (blockDim.x*gridDim.x)*threadIdx.y;
    int thread_x = threadIdx.x;
    int index = block_y + block_x + thread_y + thread_x;
    int shared_index = threadIdx.y *blockDim.x + threadIdx.x;

    // recheck this oob shit
    if (blockIdx.x*blockDim.x + threadIdx.x < rows && blockIdx.y*blockDim.y + threadIdx.y < cols) {
        s_input2[shared_index] = d_grey[index];
    } else {
        s_input2[shared_index] = 0;
    }

    __syncthreads();

    for (int r=-1; r<2; r++) {
        for (int c=-1; c<2; c++)
	{
	    int o_val = 0;
	    int sobel_y = r*(gridDim.x*blockDim.x);
	   
	   /* 
	    // if less/more than cached
	    if ((int)threadIdx.x + c < 0    || (int)threadIdx.y + r < 0 || (int)threadIdx.x + c > 15 || (int)threadIdx.y + r > 15) {
	        // if less than image
		if (block_y + thread_y + sobel_y < 0    || block_x + thread_x + c < 0 || block_y + thread_y + sobel_y > cols || block_x + thread_x + c > rows) {
		    o_val = 0;
		} else {
		// if still in bounds of image, grab val from image
		   int aug_index = index + r*(blockDim.x*gridDim.x) + c;
		    o_val =  d_grey[aug_index];
		}	
            } else {
	    // in bounds of cache, grab from cache
	    	o_val = s_input2[shared_index + sobel_y + c];
	    }
	*/
	    o_val = s_input2[shared_index + sobel_y + c];
	    sum_x = sum_x + o_val * dx[1+r][1+c];
	    sum_y = sum_y + o_val * dy[1+r][1+c];
	}
    }

    // Set your output to = (uchar)(abs(sum_x)+abs(sum_y));
    d_sobel[index] = (unsigned char)(abs(sum_x)+abs(sum_y));
}

void
your_rgba_to_greyscale(uchar4* d_rgba, unsigned char* d_grey, int rows, int cols) {
    // going to somewhat arbitrarily use 128 threads per block
    // only doing one calculation per pixel, keeping 1 dimensional
    // pixels represennted as 1d array
    // use formula (N + M-1) / M to calc grid size

    // **********************************
    // [    ][    ][    ][     ][    ][    ]

    // there will be indexing issues, account for in kernel

    int N = rows*cols;
    int M = 128;
    const dim3 blockSize(M);
    const dim3 gridSize( (N + M-1)/M );
    rgba_to_greyscale<<<gridSize,blockSize>>>(d_rgba, d_grey, N);
    cudaDeviceSynchronize();
}

void
your_sobel(unsigned char* d_grey, unsigned char* d_sobel, int rows, int cols)
{
    int N = rows*cols;
    int M = 16;
    const dim3 blockDimHist(M,M,1);
    const dim3 gridDimHist( (rows + M-1)/M, (cols + M-1)/M );
    size_t blockSharedMemory = blockDimHist.x*blockDimHist.y*sizeof(uchar);
    sobel_filter<<<gridDimHist, blockDimHist, blockSharedMemory>>>(d_grey, d_sobel, rows, cols);
    cudaDeviceSynchronize(); // maybe delete this
}

int main(int argc, char **argv)
{
    uchar4        *h_rgbaImage, *d_rgbaImage;
    unsigned char *h_greyImage, *d_greyImage;
    unsigned char *h_sobel, *d_sobel;
    std::string input_file;
    std::string greyscale_file;
    std::string sobel_file;

    switch(argc)
    {
        
	case 2:
            input_file = std::string(argv[1]);
            greyscale_file = "project1G_greyscale.png";
            sobel_file = "project1G_sobel.png";
            break;
        case 3:
            input_file = std::string(argv[1]);
            greyscale_file = std::string(argv[2]);
            sobel_file = "project1G_sobel.png";
            break;
        case 4:
            input_file = std::string(argv[1]);
            greyscale_file = std::string(argv[2]);
            sobel_file = std::string(argv[3]);
            break;
        default:
            std::cerr << "Usage: ./project1G input_file [greyscale_filename]" 
                      << "[sobel_filename]"
                      << std::endl;
            exit(1);
    }

    cv::Mat image;
    image = cv::imread(input_file.c_str(), CV_LOAD_IMAGE_COLOR);
    if(image.empty())
    {
        std::cerr << "Couldn't open file: " << input_file << std::endl;
        exit(1);
    }   
    cv::Mat imageRGBA;
    cv::Mat imageGrey;
    cv::Mat imageSobel;
    cv::cvtColor(image, imageRGBA, CV_BGR2RGBA);

    imageGrey.create(image.rows, image.cols, CV_8UC1);
    imageSobel.create(image.rows, image.cols, CV_8UC1);

    if(!imageRGBA.isContinuous() || !imageGrey.isContinuous())
    {
        std::cerr << "Images aren't continuous. Exiting" << std::endl;
        exit(1);
    }

    *(&h_rgbaImage) = (uchar4*)imageRGBA.ptr<unsigned char>(0);
    *(&h_greyImage) = imageGrey.ptr<unsigned char>(0);
    *(&h_sobel) = imageSobel.ptr<unsigned char>(0);
    size_t numRows = imageRGBA.rows;
    size_t numCols = imageRGBA.cols;
    size_t numPixels = numRows * numCols;

    // Allocate all your memory here. Use cudaMallocs, Memset, and Memcpy
    //std::cout << "Input Image: R x C: " << numRows << " x " << numCols << std::endl;

    // print rgba image pixel vals
    //std::cout << "RGBA Value Matrix:" << std::endl;
    /*
    for (int i = 1; i <= numPixels; i++) {

	    std::cout << "(" << static_cast<unsigned>(h_rgbaImage[i-1].x) << " "
		      << static_cast<unsigned>(h_rgbaImage[i-1].y) << " "
		      << static_cast<unsigned>(h_rgbaImage[i-1].z) << " "
		      << static_cast<unsigned>(h_rgbaImage[i-1].w) << ") ";
    	    if (i % numCols == 0) {
	    	std::cout << std::endl;
	    }
    }
    */

	
    // device cudaMalloc for rgba and grey images
    cudaMalloc( (uchar4 **)&d_rgbaImage, sizeof(uchar4)*numPixels );
    cudaMalloc( (unsigned char **)&d_greyImage, sizeof(unsigned char)*numPixels );

    // host to device cudaMemcpy for rgba image
    cudaMemcpy(d_rgbaImage, h_rgbaImage, sizeof(uchar4)*numPixels, cudaMemcpyHostToDevice);

    your_rgba_to_greyscale(d_rgbaImage, d_greyImage, numRows, numCols);
    cudaDeviceSynchronize();

    // Do a memcpy for your grey Image here (Device to Host)
    cudaMemcpy(h_greyImage, d_greyImage, sizeof(unsigned char)*numPixels, cudaMemcpyDeviceToHost);

    // print rgba image pixel vals
    /*
    std::cout << "GREY Value Matrix:" << std::endl;
    for (int i = 1; i <= numPixels; i++) {

	    std::cout << "(" << static_cast<unsigned>(h_greyImage[i-1]) << ") ";
    	    if (i % numCols == 0) {
	    	std::cout << std::endl;
	    }
    }
    */

    // Write out greyscale
    cv::Mat output(numRows, numCols, CV_8UC1, (void*)h_greyImage);
    cv::imwrite(greyscale_file.c_str(), output);

    cudaDeviceSynchronize();

    // Now we'll do the Sobel filter

    // cudaMalloc for sobel output
    cudaMalloc( (unsigned char**)&d_sobel, sizeof(unsigned char)*numPixels);

    your_sobel(d_greyImage, d_sobel, numRows, numCols);

    // device to host memcpy of sobel
    cudaMemcpy(h_sobel, d_sobel, sizeof(unsigned char)*numPixels, cudaMemcpyDeviceToHost);

    // Write Sobel
    cv::Mat output2(numRows, numCols, CV_8UC1, (void*)h_sobel);
    cv::imwrite(sobel_file.c_str(), output2);

    // Free \o/
    cudaFree(d_rgbaImage);
    cudaFree(d_greyImage);
    cudaFree(d_sobel);
    return 0;
}
