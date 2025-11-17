#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "bvh.h"
#include "ray.h"
#include "camera.h"
#include "pointlight.h"
#include "hit_struct.h"
#include <string>
#include <vector>

using namespace std;

class Raytracer{
    public:

        static vector<float> reflectRefract(Ray, BVHNode*, const Camera&, const vector<PointLight>&, int = 0);
        static vector<float> blinnPhong(const HitStructure&, const Camera&, const vector<PointLight>);

        //Helper functions
        static vector<float> normalise(vector<float>);
        static float dotProd(vector<float>, vector<float>);
        static vector<float> add_vec(vector<float>, vector<float>);
        static vector<float> sub_vec(vector<float>, vector<float>);
        static vector<float> mul_vec(vector<float>, float);
};
#endif