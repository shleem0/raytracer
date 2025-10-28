#ifndef PLANE_H
#define PLANE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include <vector>
#include <string>

using namespace std;

class Plane : public Shape{
    public:
        vector<vector<float>> vertices;
        float normal[3];

        static vector<Plane> parsePlaneDataFromJson();
        bool intersect(const Ray&, HitStructure&) override;

        bool pointInsidePolygon(const vector<float>&);
        void calculateNormal();
        void sortVerticesWinding();

        AABB getAABB() const;
};

#endif