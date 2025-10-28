#ifndef SHAPE
#define SHAPE

#include "ray.h"
#include "hit_struct.h"

using namespace std;

class Shape{
    public:
        float location[3];
        static float getFloat(const string&, const string&);
        static string getJSONObject(const string&, const string&);
        static string getJSONArray(const string&, const string&);
        virtual bool intersect(const Ray&, HitStructure&) = 0;
};

#endif