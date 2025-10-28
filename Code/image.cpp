#include <iostream>
#include <stdexcept>
#include <string>
#include "image.h"
using namespace std;

Image::Image(const string& fileName){
    file = fileName;
}

//Function to read image, pixel by pixel
void Image::readImage() {
    ifstream fileStream("../Textures/" + file);

    if (!fileStream.is_open()) {
        throw runtime_error("Failed to open the image file");
    }

    string header;
    fileStream >> header;

    if (header != "P3" && header != "P6") {
        throw runtime_error("Only P3 and P6 .ppm files are supported");
    }

    fileStream >> width >> height >> maxColourValue;

    pixelValues.resize(height, vector<vector<int>>(width, vector<int>(3)));

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

    fileStream.close();
}

//Function to modify pixels colour by RGB values
void Image::modifyPixel(int x, int y, int red, int green, int blue) {
    if (x < 0 || x >= width || y < 0 || y >= height) {
        throw runtime_error("Pixel coordinates are out of bounds");
    }

    pixelValues[y][x][0] = red;
    pixelValues[y][x][1] = green;
    pixelValues[y][x][2] = blue;
}


//Function to write image to chosen file
void Image::writeImage(const string& outputFileName) {
    ofstream fileStream("../Textures/" + outputFileName);

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