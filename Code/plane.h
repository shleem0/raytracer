#ifndef PLANE_H
#define PLANE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include "config.h"
#include <vector>
#include <string>

using namespace std;

//Class to represent a plane in a Blender scene
class Plane : public Shape{
    public:
        vector<vector<float>> vertices;
        vector<float> normal;

        //Parse all plane data from the scene's JSON
        static vector<Plane> parsePlaneDataFromJson(Config);

        //Check if a ray intersects with the plane
        bool intersect(const Ray&, HitStructure&, Config) override;

        //Check if a 2D point is within the plane's bounds (as a 2D polygon)
        bool pointInsidePolygon(const vector<float>&);

        //Calculate normal of the plane
        void calculateNormal();

        //Sort the planes by angle around a centroid
        void sortVerticesWinding();

        //Get AABB of plane
        AABB getAABB() const;
};

#endif