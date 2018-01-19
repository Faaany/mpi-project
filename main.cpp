#include <iostream>
#include <vector>
#include "lodepng.h"
#include <string>
#include <mpi.h>

using std::string;

typedef std::vector<double> Vector;

Vector LoadImage(const string& filename, unsigned& width, unsigned& height);
void SaveImage(const Vector &vector, unsigned width, unsigned height, const string &filename);

int main() {

    unsigned width, height;

    Vector image = LoadImage("image.png", width, height);

    std::cout << "The dimension of the image is " << height << "x" << width << ".\n";

    int id, num_procs;

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    


    MPI_Finalize();



    SaveImage(image, width, height,"result.png");

    return 0;
}

void SaveImage(const Vector &image, unsigned width, unsigned height, const string &filename) {
    std::vector<unsigned char> pixels(width*height*4);
    for(unsigned i=0; i<width*height; ++i) {
        pixels[4*i] = pixels[4*i+1] = pixels[4*i+2]=(unsigned char)image[i];
        pixels[4*i+3]= 255;
    }

    unsigned error = lodepng::encode(filename, pixels, width, height);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

}


Vector LoadImage(const string& filename, unsigned& width, unsigned& height) {
    std::vector<unsigned char> image; //the raw pixels
    width=0, height=0;

    //decode
    unsigned error = lodepng::decode(image, width, height, filename);

    //if there's an error, display it
    if(error) std::cerr << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;

    Vector result(width*height);
    for(unsigned i=0; i<width*height; ++i) {
        result[i]= 1.0/3*(image[4*i]+image[4*i+1]+image[4*i+2]);
    }

    return result;
}
