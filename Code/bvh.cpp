#include <vector>
#include <iostream>
#include "bvh.h"

using namespace std;

//Build the BVH node for a given AABB
BVHNode::BVHNode(const AABB& box) : aabb(box), left(nullptr), right(nullptr) {}

//Build the next BVH node given a child BVH node
BVHNode::BVHNode(const BVHNode& other)
    : aabb(other.aabb), left(nullptr), right(nullptr) {
    planes = other.planes;
    cubes = other.cubes;
    spheres = other.spheres;
    if (other.left) {
        left = new BVHNode(*other.left);
    }
    if (other.right) {
        right = new BVHNode(*other.right);
    }
}

BVHNode::~BVHNode() {
    if (left) {
        delete left;
    }
    if (right) {
        delete right;
    }
}


bool BVHNode::intersect(Ray& ray, float& tMin, float& tMax){
    //Check if ray intersects this node's bounding box
    if (!aabb.intersect(ray, tMin, tMax)) {
        return false;
    }

    bool hit = false;
    float tLeftMin = tMin, tLeftMax = tMax;
    float tRightMin = tMin, tRightMax = tMax;

    //Recursively check left child
    if (left && left->intersect(ray, tLeftMin, tLeftMax)) {
        tMin = tLeftMin;
        tMax = tLeftMax;
        hit = true;
    }

    //Recursively check right child
    if (right && right->intersect(ray, tRightMin, tRightMax)) {
        if (!hit || tRightMin < tMin) {
            tMin = tRightMin;
            tMax = tRightMax;
        }
        hit = true;
    }

    //If hit an internal child, return true
    if (hit) return true;

    //Leaf node - check primitives
    if (!left && !right) {
        for (auto plane : planes) {
            HitStructure hs;
            if (plane.intersect(ray, hs)) {
                hs.objectType = "Plane";
                ray.hs.push_back(hs);
                return true;
            }
        }

        for (auto cube : cubes) {
            HitStructure hs;
            if (cube.intersect(ray, hs)) {
                hs.objectType = "Cube";
                ray.hs.push_back(hs);
                return true;
            }
        }

        for (auto sphere : spheres) {
            HitStructure hs;
            if (sphere.intersect(ray, hs)) {
                hs.objectType = "Sphere";
                ray.hs.push_back(hs);
                return true;
            }
        }
    }

    return false; //No intersection
}


BVHNode* BVHNode::buildBVH(const vector<Plane> planes, const vector<Cube> cubes, const vector<Sphere> spheres, int maxDepth) {
    if (maxDepth == 0) {
        //Create a leaf node
        AABB aabb = AABB::fromPoints(planes, cubes, spheres);
        BVHNode* node = new BVHNode(aabb);
        node->planes.insert(node->planes.end(), planes.begin(), planes.end());
        node->cubes.insert(node->cubes.end(), cubes.begin(), cubes.end());
        node->spheres.insert(node->spheres.end(), spheres.begin(), spheres.end());
        return node;
    }

    //Calculate the bounding box for all objects
    AABB aabb = AABB::fromPoints(planes, cubes, spheres);

    //Choose the axis to split on (x, y, or z)
    int axis = 0;  //X-axis
    if (aabb.max[1] - aabb.min[1] > aabb.max[0] - aabb.min[0]) {
        axis = 1;  //Y-axis
    }
    if (aabb.max[2] - aabb.min[2] > aabb.max[0] - aabb.min[0] && aabb.max[2] - aabb.min[2] > aabb.max[1] - aabb.min[1]) {
        axis = 2;  //Z-axis
    }

    //Split the objects into two sets
    vector<Plane> leftPlanes;
    vector<Cube> leftCubes;
    vector<Sphere> leftSpheres;
    vector<Plane> rightPlanes;
    vector<Cube> rightCubes;
    vector<Sphere> rightSpheres;

    float splitPoint = (aabb.min[axis] + aabb.max[axis]) / 2.0f;
    for (const auto& plane : planes) {
        if ((plane.getAABB().min[axis] + plane.getAABB().max[axis]) / 2.0f < splitPoint) {
            leftPlanes.push_back(plane);
        } else {
            rightPlanes.push_back(plane);
        }
    }
    for (const auto& cube : cubes) {
        if ((cube.getAABB().min[axis] + cube.getAABB().max[axis]) / 2.0f < splitPoint) {
            leftCubes.push_back(cube);
        } else {
            rightCubes.push_back(cube);
        }
    }
    for (const auto& sphere : spheres) {
        if ((sphere.getAABB().min[axis] + sphere.getAABB().max[axis]) / 2.0f < splitPoint) {
            leftSpheres.push_back(sphere);
        } else {
            rightSpheres.push_back(sphere);
        }
    }

    //Recursively build the left and right child nodes
    BVHNode* node = new BVHNode(aabb);
    if (!leftPlanes.empty() || !leftCubes.empty() || !leftSpheres.empty()) {
        node->left = buildBVH(leftPlanes, leftCubes, leftSpheres, maxDepth - 1);
    }
    if (!rightPlanes.empty() || !rightCubes.empty() || !rightSpheres.empty()) {
        node->right = buildBVH(rightPlanes, rightCubes, rightSpheres, maxDepth - 1);
    }

    return node;
}
