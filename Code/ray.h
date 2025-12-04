#ifndef RAY
#define RAY

#include "hit_struct.h"
#include <vector>

//Class to represent a ray in 3D space
class Ray{
    public:
        vector<float> origin;
        vector<float> direction;
        vector<HitStructure> hs;
        float time = 0.0f;
};
#endif