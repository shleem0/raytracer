#ifndef IMAGE
#define IMAGE

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

//Class representing a P3 or P6 image
class Image{
    public:
        string file;
        int width;
        int height;
        int maxColourValue;
        vector<vector<vector<int>>> pixelValues;

        //Constructor taking file name
        Image(const string&);

        //Read the image pixel by pixel from the file
        void readImage();

        //Change a pixel by giving new RGB values for its colour
        void modifyPixel(int, int, int, int, int);

        //Write the image to a ppm file
        void writeImage(const string&);
};

#endif