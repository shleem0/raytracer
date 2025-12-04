#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "bvh.h"
#include "ray.h"
#include "camera.h"
#include "pointlight.h"
#include "hit_struct.h"
#include "config.h"
#include <string>
#include <vector>

using namespace std;

class Raytracer{
    public:

        //Main functions
        static vector<float> reflectRefract(Config, BVHNode*, Ray&, const Camera&, vector<PointLight>&, vector<Plane>&, vector<Cube>&, vector<Sphere>&, int = 0, int = 5);
        static vector<float> blinnPhong(Config, BVHNode*, HitStructure&, const Camera&, vector<PointLight>&, vector<Plane>&, vector<Cube>&, vector<Sphere>&);
        static float computeHardShadows(const vector<float>&, const float, const vector<float>&, const PointLight&, BVHNode*, vector<Plane>&, vector<Cube>&, vector<Sphere>&, Config);
        static bool unacceleratedIntersection(Ray&, vector<Plane>&, vector<Sphere>&, vector<Cube>&, Config);

        //Distributed RT
        static float computeSoftShadows(const vector<float>&, const float, const vector<float>&, const PointLight&, BVHNode*, vector<Plane>&, vector<Cube>&, vector<Sphere>&, Config);

        //Helper functions
        //Vectors
        static vector<float> normalise(vector<float>);
        static float dotProd(vector<float>, vector<float>);
        static vector<float> crossProd(vector<float>, vector<float>);
        static vector<float> add_vec(vector<float>, vector<float>);
        static vector<float> sub_vec(vector<float>, vector<float>);
        static vector<float> mul_vec(vector<float>, float);

        //Random sampling
        static vector<float> rndUnitSphere();
        static vector<float> rndConeDirection(const vector<float>&, float);

        //Other
        static float lerp(float, float, float);
};
#endif