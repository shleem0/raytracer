#ifndef CUBE_H
#define CUBE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include <string>

using namespace std;

//Class to represent a cube in a Blender scene
class Cube : public Shape{
    public:
        vector<float> rotation;
        float scale;

        //Method to extract all cube data from the scene's JSON
        static vector<Cube> parseCubeDataFromJson();

        //Checks if a ray intersects with a cube
        bool intersect(const Ray&, HitStructure&) override;

        //Rotation transformations for intersection checking
        void rotateXYZInverse(vector<float>&, vector<float>&);
        void rotateXYZ(vector<float>&, vector<float>&);

        AABB getAABB() const;
};

#endif