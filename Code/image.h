#ifndef IMAGE
#define IMAGE

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

class Image{
    public:
        string file;
        int width;
        int height;
        int maxColourValue;
        vector<vector<vector<int>>> pixelValues;

        Image(const string&);
        void readImage();
        void modifyPixel(int, int, int, int, int);
        void writeImage(const string&);
};

#endif