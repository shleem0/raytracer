#ifndef RAY
#define RAY

#include "hit_struct.h"
#include <vector>

class Ray{
    public:
        float origin[3];
        float direction[3];
        vector<HitStructure> hs;
};
#endif