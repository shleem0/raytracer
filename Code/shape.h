#ifndef SHAPE
#define SHAPE

#include "ray.h"
#include "hit_struct.h"
#include "image.h"
#include <vector>
#include <string>

using namespace std;

//Superclass to represent a generic shape in a Blender scene
class Shape{
    public:
        vector<float> location;
        vector<float> diffuse;
        vector<float> specular;
        float shininess;
        float transparency;
        float ior;
        string texture;
        bool hasTex = false;

        //Parse a float value from a JSON object
        static float getFloat(const string&, const string&);

        static string getString(const string&, const string&);

        //Parse a string value from a JSON object
        static string getJSONObject(const string&, const string&);

        //Parse an array from a JSON object
        static string getJSONArray(const string&, const string&);

        //Abstract ray intersection class implemented by each shape
        virtual bool intersect(const Ray&, HitStructure&) = 0;
};

#endif