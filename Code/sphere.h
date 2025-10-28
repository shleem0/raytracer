#ifndef SPHERE_H
#define SPHERE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include <string>

using namespace std;

class Sphere : public Shape{
    public:
        float radius;

        static vector<Sphere> parseSphereDataFromJson();
        bool intersect(const Ray&, HitStructure&) override;

        AABB getAABB() const;
};

#endif