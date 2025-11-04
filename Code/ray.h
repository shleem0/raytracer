#ifndef RAY
#define RAY

#include "hit_struct.h"
#include <vector>

//Class to represent a ray in 3D space
class Ray{
    public:
        float origin[3];
        float direction[3];
        vector<HitStructure> hs;
};
#endif