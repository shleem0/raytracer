#include "camera.h"
#include "image.h"
#include "ray.h"
#include "shape.h"
#include "sphere.h"
#include "cube.h"
#include "plane.h"
#include "hit_struct.h"
#include "aabb.h"
#include "bvh.h"
#include <iostream>
#include <ctime>
#include <vector>
#include <chrono>
#include <cfloat>
using namespace std;
using namespace std::chrono;

int main(){

    vector<Camera> cams = Camera::parseCameraDataFromJson();
    vector<Sphere> spheres = Sphere::parseSphereDataFromJson();
    vector<Cube> cubes = Cube::parseCubeDataFromJson();
    vector<Plane> planes = Plane::parsePlaneDataFromJson();

    int totalObjects = spheres.size() + cubes.size() + planes.size();

    bool intersect;
    Ray ray;
    vector<Ray> rays;

    cout << "---Testing for " << totalObjects <<" objects---\n";


    /*Accelerated ray caster*/    
    //Build BVH
    BVHNode* bvh = BVHNode::buildBVH(planes, cubes, spheres);

    uint64_t start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    //For every pixel
    for (int x = 0; x < cams[0].res_x; x++) {
        for (int y = 0; y < cams[0].res_y; y++) {
            intersect = false;

            ray = cams[0].convertPixelToRay(x, y);
            rays.push_back(ray);

            //Traverse BVH
            float tMin = 0.0f;
            float tMax = FLT_MAX;
            if (bvh->intersect(ray, tMin, tMax)) {
                intersect = true;
            }
        }
    }
    uint64_t end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

    cout << "Time taken for BVH RT (accelerated): " << (end - start) / 1000.0f << "s\n";

    delete bvh;

    //---Simple ray caster (unaccelerated)---
    start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

    //For every pixel
    for (int x = 0; x < cams[0].res_x; x++){
        for (int y = 0; y < cams[0].res_y; y++){

            intersect = false;

            ray = cams[0].convertPixelToRay(x, y);
            rays.push_back(ray);

            //Check intersection with planes
            for (int p = 0; p < planes.size(); p++){

                HitStructure planeHS;
                intersect = planes[p].intersect(ray, planeHS);

                if (intersect){
                    ray.hs.push_back(planeHS);
                    planeHS.objectType = "Plane";
                    break;
                }
            }

            //Check intersection with cubes
            for (int c = 0; c < cubes.size(); c++){

                HitStructure cubeHS;
                intersect = cubes[c].intersect(ray, cubeHS);

                if (intersect){
                    ray.hs.push_back(cubeHS);
                    cubeHS.objectType = "Cube";
                    break;
                }
            }

            //Check intersection with spheres
            for (int s = 0; s < spheres.size(); s++){

                HitStructure sphereHS;
                intersect = spheres[s].intersect(ray, sphereHS);

                if (intersect){
                    ray.hs.push_back(sphereHS);
                    sphereHS.objectType = "Sphere";
                    break;
                }
            }
        }
    }

    end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    cout << "Time taken for standard RT (non-accelerated): " << (end - start) / 1000.0f << "s\n";

    return 0;
}

