#ifndef HIT_STRUCT
#define HIT_STRUCT

#include <string>
#include <vector>

using namespace std;

//Data on a ray intersectin with a shape
struct HitStructure{
    vector<float> hitPoint;
    float rayDistance;

    vector<float> normal;
    vector<float> diffuse;
    vector<float> specular;
    float shininess;
    float transparency;
    float ior;
    string objectType;
};

#endif