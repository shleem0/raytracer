#ifndef SPHERE_H
#define SPHERE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include <string>

using namespace std;

//Class to represent a sphere in a Blender scene
class Sphere : public Shape{
    public:
        float radius;

        //Parse all sphere data from the scene's JSON
        static vector<Sphere> parseSphereDataFromJson();

        //Check if a ray intersects with a sphere
        bool intersect(const Ray&, HitStructure&) override;

        //Get sphere's AABB
        AABB getAABB() const;
};

#endif