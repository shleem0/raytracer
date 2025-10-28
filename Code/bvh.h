#ifndef BVH
#define BVH

#include <vector>
#include "ray.h"
#include "aabb.h"
#include "plane.h"
#include "cube.h"
#include "sphere.h"

class BVHNode {
public:
    AABB aabb;
    BVHNode* left;
    BVHNode* right;
    vector<Plane> planes;
    vector<Cube> cubes;
    vector<Sphere> spheres;

    // Constructor
    BVHNode(const AABB& aabb);
    BVHNode(const BVHNode& other);
    ~BVHNode();

    // Traversing the BVH for ray intersection
    bool intersect(Ray& ray, float& tMin, float& tMax);

    // Build the BVH
    static BVHNode* buildBVH(const vector<Plane> planes, const vector<Cube> cubes, const vector<Sphere> spheres, int maxDepth = 8);
};

#endif 
