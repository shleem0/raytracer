#ifndef BVH
#define BVH

#include <vector>
#include "ray.h"
#include "aabb.h"
#include "plane.h"
#include "cube.h"
#include "sphere.h"
#include "config.h"

//Class for a single node in the bounding volume hierarchy
class BVHNode {
public:
    //Current AABB
    AABB aabb;

    //Left and right children of the current node
    BVHNode* left;
    BVHNode* right;

    //Shapes in the scene
    vector<Plane> planes;
    vector<Cube> cubes;
    vector<Sphere> spheres;

    //Constructor and destructor
    BVHNode(const AABB& aabb);
    BVHNode(const BVHNode& other);
    ~BVHNode();

    //Traversing the BVH for ray intersection
    bool intersect(Ray& ray, float& tMin, float& tMax, Config config);

    //Build the BVH
    static BVHNode* buildBVH(const vector<Plane> planes, const vector<Cube> cubes, const vector<Sphere> spheres, Config config, int maxDepth = 8);
};

#endif 
