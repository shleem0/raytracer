#ifndef HIT_STRUCT
#define HIT_STRUCT

#include <string>
#include <vector>

using namespace std;

struct HitStructure{
    vector<float> hitPoint;
    float rayDistance;
    string objectType;
};

#endif