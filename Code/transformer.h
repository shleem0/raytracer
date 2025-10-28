#ifndef TRANSFORMER
#define TRANSFORMER

#include "ray.h"
#include "cube.h"

using namespace std;

class Transformer{
    public:
        static void rotateRayInverse(const Ray&, Ray&, Cube*);
        static void transRayToLocal(const Ray&, Ray&, Cube*);
        static bool intersectUnitCube(const Ray&, float, float);

        static void rotateXYZInverse(vector<float>);
        static void rotateXYZ(vector<float>);
};

#endif