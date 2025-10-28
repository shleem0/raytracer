#ifndef CUBE_H
#define CUBE_H

#include "shape.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include <string>

using namespace std;



class Cube : public Shape{
    public:
        float rotation[3];
        float scale;

        static vector<Cube> parseCubeDataFromJson();
        bool intersect(const Ray&, HitStructure&) override;

        void rotateXYZInverse(vector<float>&, float*);
        void rotateXYZ(vector<float>&, float*);

        AABB getAABB() const;
};

#endif