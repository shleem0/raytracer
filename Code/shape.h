#ifndef SHAPE
#define SHAPE

#include "ray.h"
#include "hit_struct.h"

using namespace std;

//Superclass to represent a generic shape in a Blender scene
class Shape{
    public:
        float location[3];
        float diffuse[3];
        float specular[3];
        float shininess;
        float transparency;
        float ior;

        //Parse a float value from a JSON object
        static float getFloat(const string&, const string&);

        //Parse a string value from a JSON object
        static string getJSONObject(const string&, const string&);

        //Parse an array from a JSON object
        static string getJSONArray(const string&, const string&);

        //Abstract ray intersection class implemented by each shape
        virtual bool intersect(const Ray&, HitStructure&) = 0;
};

#endif