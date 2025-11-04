#include "aabb.h"
#include "plane.h"
#include "cube.h"
#include "sphere.h"
#include "ray.h"
#include <vector>
#include <cmath>


using namespace std;

//Constructor
AABB::AABB(vector<float> min_v, vector<float> max_v){
    min = min_v;
    max = max_v;

}

AABB::AABB(const AABB& other) : min(other.min), max(other.max) {}


//Checks if a point is within the AABB
bool AABB::contains(const vector<float> point) const {
    return (point[0] >= min[0] && point[0] <= max[0]) &&
           (point[1] >= min[1] && point[1] <= max[1]) &&
           (point[2] >= min[2] && point[2] <= max[2]);
}

//Checks if the given ray intersects with the AABB
bool AABB::intersect(const Ray& ray, float& tMin, float& tMax) const {
    float tMinX, tMaxX, tMinY, tMaxY, tMinZ, tMaxZ;

    if (ray.direction[0] >= 0) {
        tMinX = (min[0]- ray.origin[0]) / ray.direction[0];
        tMaxX = (max[0] - ray.origin[0]) / ray.direction[0];
    } else {
        tMinX = (max[0] - ray.origin[0]) / ray.direction[0];
        tMaxX = (min[0] - ray.origin[0]) / ray.direction[0];
    }

    if (ray.direction[1] >= 0) {
        tMinY = (min[1] - ray.origin[1]) / ray.direction[1];
        tMaxY = (max[1] - ray.origin[1]) / ray.direction[1];
    } else {
        tMinY = (max[1] - ray.origin[1]) / ray.direction[1];
        tMaxY = (min[1] - ray.origin[1]) / ray.direction[1];
    }

    if (ray.direction[2] >= 0) {
        tMinZ = (min[2] - ray.origin[2]) / ray.direction[2];
        tMaxZ = (max[2] - ray.origin[2]) / ray.direction[2];
    } else {
        tMinZ = (max[2] - ray.origin[2]) / ray.direction[2];
        tMaxZ = (min[2] - ray.origin[2]) / ray.direction[2];
    }
    

    tMin = std::max(tMinX, std::max(tMinY, tMinZ));
    tMax = std::min(tMaxX, std::min(tMaxY, tMaxZ));

    return tMin <= tMax;
}

//Gets AABBs of each shape in the scene
AABB AABB::fromPoints(const vector<Plane> planes, const vector<Cube> cubes, const vector<Sphere> spheres) {
    vector<float> minPoint = {numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()};
    vector<float> maxPoint = {numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min()};

    for (const auto& plane : planes) {
        AABB aabb = plane.getAABB();
        for (int i = 0; i < 3; i++) {
            if (aabb.min[i] < minPoint[i]) {
                minPoint[i] = aabb.min[i];
            }
            if (aabb.max[i] > maxPoint[i]) {
                maxPoint[i] = aabb.max[i];
            }
        }
    }

    for (const auto& cube : cubes) {
        AABB aabb = cube.getAABB();
        for (int i = 0; i < 3; i++) {
            if (aabb.min[i] < minPoint[i]) {
                minPoint[i] = aabb.min[i];
            }
            if (aabb.max[i] > maxPoint[i]) {
                maxPoint[i] = aabb.max[i];
            }
        }
    }

    for (const auto& sphere : spheres) {
        AABB aabb = sphere.getAABB();
        for (int i = 0; i < 3; i++) {
            if (aabb.min[i] < minPoint[i]) {
                minPoint[i] = aabb.min[i];
            }
            if (aabb.max[i] > maxPoint[i]) {
                maxPoint[i] = aabb.max[i];
            }
        }
    }

    return AABB(minPoint, maxPoint);
}


//Get AABB from sent of points around a shape
AABB AABB::fromPoints(const vector<vector<float>> points) {
    vector<float> minPoint = {numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()};
    vector<float> maxPoint = {numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min()};

    for (const auto& point : points) {
        for (int i = 0; i < 3; i++) {
            if (point[i] < minPoint[i]) {
                minPoint[i] = point[i];
            }
            if (point[i] > maxPoint[i]) {
                maxPoint[i] = point[i];
            }
        }
    }

    return AABB(minPoint, maxPoint);
}


//Get the union of two AABBs
AABB AABB::unionOf(const AABB& a, const AABB& b) {

    vector<float> minPoint = {numeric_limits<float>::max(), numeric_limits<float>::max(), numeric_limits<float>::max()};
    vector<float> maxPoint = {numeric_limits<float>::min(), numeric_limits<float>::min(), numeric_limits<float>::min()};

    for (int i = 0; i < 3; i++) {
        if (a.min[i] < minPoint[i]) {
            minPoint[i] = a.min[i];
        }
        if (b.min[i] < minPoint[i]) {
            minPoint[i] = b.min[i];
        }
        if (a.max[i] > maxPoint[i]) {
            maxPoint[i] = a.max[i];
        }
        if (b.max[i] > maxPoint[i]) {
            maxPoint[i] = b.max[i];
        }
    }

    return AABB(minPoint, maxPoint);
}
