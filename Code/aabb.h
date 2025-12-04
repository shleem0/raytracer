#ifndef AABB_H
#define AABB_H

#include <vector>
#include "config.h"

using namespace std;

class Plane;
class Cube;
class Sphere;
class Ray;


//The class representing axis-aligned bounding boxes
class AABB{
    public:
        vector<float> min;  //Minimum point (lower-left-front corner)
        vector<float> max;  //Maximum point (upper-right-back corner)

        //Constructor
        AABB(vector<float> minPoint, vector<float> maxPoint);
        AABB(const AABB& other);

        //Check if a point is inside the AABB
        bool contains(const vector<float>) const;

        //Check if a ray intersects the AABB
        bool intersect(const Ray&, float&, float&) const;

        //Calculate the AABB for all shapes in the scene
        static AABB fromPoints(const vector<Plane>, const vector<Cube>, const vector<Sphere>, Config);

        //Calculate the AABB for a list of points
        static AABB fromPoints(const vector<vector<float>>);

        //Union of two AABBs
        static AABB unionOf(const AABB&, const AABB&);
};

#endif
