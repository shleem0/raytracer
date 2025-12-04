#include <iostream>
#include <stdexcept>
#include <string>
#include "image.h"
using namespace std;


Image::Image(const string& fileName){
    file = fileName;
}


void Image::readImage() {

    ifstream fileStream("../Textures/" + file, ios::binary);

    if (!fileStream.is_open()) {
        throw runtime_error("Failed to open the image file");
    }

    string header;
    fileStream >> header;

    if (header != "P3" && header != "P6") {
        throw runtime_error("Only P3 and P6 .ppm files are supported");
    }

    fileStream >> width >> height >> maxColourValue;

    //Adjust pixel colour array to size of image
    pixelValues.resize(height, vector<vector<int>>(width, vector<int>(3)));

    if (header == "P3"){
        //ASCII PPM
        //Store each pixel's colour
        int value;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                for (int color = 0; color < 3; ++color) {
                    fileStream >> value;
                    if (fileStream.fail()) {
                        throw runtime_error("Invalid pixel values");
                    }
                    pixelValues[y][x][color] = value;
                }
            }
        }
    }

    else{
        // Binary PPM
        char r, g, b;
        for (int y = 0; y < height; ++y){
            for (int x = 0; x < width; ++x) {
                fileStream.read(&r, 1);
                fileStream.read(&g, 1);
                fileStream.read(&b, 1);
                if (fileStream.fail()) throw runtime_error("Invalid pixel values");
                pixelValues[y][x][0] = static_cast<unsigned char>(r);
                pixelValues[y][x][1] = static_cast<unsigned char>(g);
                pixelValues[y][x][2] = static_cast<unsigned char>(b);
            }
        }
    }

    fileStream.close();
}

//Adjust specific pixel colour
void Image::modifyPixel(int x, int y, int red, int green, int blue) {

    if (x < 0 || x >= width || y < 0 || y >= height) {
        throw runtime_error("Pixel coordinates are out of bounds");
    }

    pixelValues[y][x][0] = red;
    pixelValues[y][x][1] = green;
    pixelValues[y][x][2] = blue;
}


//Write image to file
void Image::writeImage(const string& outputFileName) {
    ofstream fileStream("../Output/" + outputFileName);

    if (!fileStream.is_open()) {
        throw runtime_error("Failed to open the output file");
    }

    fileStream << "P3\n" << width << " " << height << "\n" << maxColourValue << "\n";

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fileStream << pixelValues[y][x][0] << " " << pixelValues[y][x][1] << " " << pixelValues[y][x][2] << " ";
        }
        fileStream << "\n";
    }
    fileStream.close();
}