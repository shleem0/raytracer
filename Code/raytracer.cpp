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
#include "pointlight.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <chrono>
#include <cfloat>
using namespace std;
using namespace std::chrono;


//Blinn-Phong shading for a given ray at a pixel
vector<float> blinnPhong(const HitStructure& hs, const Camera& cam, const vector<PointLight> lights){

    float ka = 0.5f; //Ambient reflect
    float kd = 0.2126f * hs.diffuse[0] + 0.7152f * hs.diffuse[1] + 0.0722f * hs.diffuse[2]; //Diffuse reflect
    float ks = 0.2126f * hs.specular[0] + 0.7152f * hs.specular[1] + 0.0722f * hs.specular[2]; //Specular reflect
    float shiny = hs.shininess;

    vector<float> matColour = hs.diffuse;
    vector<float> N = hs.normal;

    vector<float> V = {cam.location[0] - hs.hitPoint[0], cam.location[1] - hs.hitPoint[1], cam.location[2] - hs.hitPoint[2]};
    float V_len = sqrt(pow(V[0], 2) + pow(V[1], 2) + pow(V[2], 2));

    V[0] /= V_len;
    V[1] /= V_len;
    V[2] /= V_len;

    vector<float> colour = {ka * matColour[0], ka * matColour[1], ka * matColour[2]};

    for (const auto& light : lights){

        vector<float> L = {light.location[0] - hs.hitPoint[0], light.location[1] - hs.hitPoint[1], light.location[2] - hs.hitPoint[2]};
        float L_len = sqrt(pow(L[0], 2) + pow(L[1], 2) + pow(L[2], 2));

        L[0] /= L_len;
        L[1] /= L_len;
        L[2] /= L_len;

        float irr = light.rad_intensity / (1.0f + pow(L_len, 2));

        vector<float> H = {L[0] + V[0], L[1] + V[1], L[2] + V[2]};
        float H_len = sqrt(pow(H[0], 2) + pow(H[1], 2) + pow(H[2], 2));

        H[0] /= H_len;
        H[1] /= H_len;
        H[2] /= H_len;

        float N_L = max(0.0f, N[0]*L[0] + N[1]*L[1] + N[2]*L[2]);
        float N_H = max(0.0f, N[0]*H[0] + N[1]*H[1] + N[2]*H[2]);

        //float diffuse = kd * N_L * irr;
        //float specular = 0.5f * ks * pow(N_H, shiny) * irr;

        colour[0] += hs.diffuse[0] * N_L * irr + 0.5f * hs.specular[0] * pow(N_H, shiny) * irr;
        colour[1] += hs.diffuse[1] * N_L * irr + 0.5f * hs.specular[1] * pow(N_H, shiny) * irr;
        colour[2] += hs.diffuse[2] * N_L * irr + 0.5f * hs.specular[2] * pow(N_H, shiny) * irr;
    }

    colour[0] = min(1.0f, max(colour[0], 0.0f));
    colour[1] = min(1.0f, max(colour[1], 0.0f));
    colour[2] = min(1.0f, max(colour[2], 0.0f));

    return colour;
}


vector<float> reflectRefract(Ray ray, BVHNode* bvh, const Camera& cam, const vector<PointLight>& lights, int depth = 0){

    const int MAX_DEPTH = 5;
    const float EPS = 1e-4f;
    const vector<float> background = {0.0f, 0.0f, 0.0f};
    float tMin = 0.0f, tMax = FLT_MAX;

    if (depth > MAX_DEPTH || !bvh->intersect(ray, tMin, tMax)){
        return background;
    }

    HitStructure hs = ray.hs.back();

    vector<float> colour = blinnPhong(hs, cam, lights);

    float reflectivity = ((hs.specular[0] + hs.specular[1] + hs.specular[2]) / 3.0f) * 0.5;
    reflectivity = min(max(reflectivity, 0.0f), 1.0f);

    float ior = hs.ior > 0.0f ? hs.ior : 1.0f;
    float transparency = hs.transparency;

    vector<float> V = {cam.location[0] - hs.hitPoint[0], cam.location[1] - hs.hitPoint[1], cam.location[2] - hs.hitPoint[2]};
    float V_len = sqrt(pow(V[0], 2) + pow(V[1], 2) + pow(V[2], 2));
    V[0] /= V_len;
    V[1] /= V_len;
    V[2] /= V_len;

    vector<float> N = hs.normal;

    float N_rayDir = ray.direction[0]*N[0] + ray.direction[1]*N[1] + ray.direction[2]*N[2];
    if (N_rayDir > 0.0f) {
        // flip normal so N faces opposite to ray.direction
        N = {-N[0], -N[1], -N[2]};
    }


    //Reflection
    vector<float> reflectColour = {0.0f, 0.0f, 0.0f};
    if (reflectivity > 0){

        float incoming_N = ray.direction[0]*N[0] + ray.direction[1]*N[1] + ray.direction[2]*N[2];
    

        vector<float> reflectDir = {ray.direction[0] - 2.0f * incoming_N * N[0],
        ray.direction[1] - 2.0f * incoming_N * N[1],
        ray.direction[2] - 2.0f * incoming_N * N[2]
        };

        float rd_len = sqrt(pow(reflectDir[0], 2) + pow(reflectDir[1], 2) + pow(reflectDir[2], 2));
        reflectDir[0] /= rd_len;
        reflectDir[1] /= rd_len;
        reflectDir[2] /= rd_len;


        vector<float> reflectOrigin = {hs.hitPoint[0] + N[0] * EPS,
        hs.hitPoint[1] + N[1] * EPS,
        hs.hitPoint[2] + N[2] * EPS};

        Ray reflectRay;
        reflectRay.origin = reflectOrigin;
        reflectRay.direction = reflectDir;
        reflectRay.hs.clear();

        float rtMin = 0.0f, rtMax = FLT_MAX;
        if (bvh->intersect(reflectRay, rtMin, rtMax)) {
            reflectColour = reflectRefract(reflectRay, bvh, cam, lights, depth + 1);
        } 
        else {
            reflectColour = background;
        }
    }

    colour[0] = colour[0] + (reflectColour[0] * reflectivity);
    colour[1] = colour[1] + (reflectColour[1] * reflectivity);    
    colour[2] = colour[2] + (reflectColour[2] * reflectivity);
    

    //Refraction
    bool refracted = false;
    if (transparency > 0){
        float n1 = 1.0f;
        float n2 = ior;
        float cos_i = ray.direction[0]*N[0] + ray.direction[1]*N[1] + ray.direction[2]*N[2];

        if (cos_i > 0){
            //N = {-N[0], -N[1], -N[2]};
            swap(n1, n2);
            //cos_i = ray.direction[0]*N[0] + ray.direction[1]*N[1] + ray.direction[2]*N[2];
        }

        float eta = n1/n2;
        cos_i = -max(-1.0f, min(1.0f, cos_i));
        vector<float> refractDir;

        float k = 1.0f - eta * eta * (1.0f - cos_i * cos_i);
        if (k < 0.0f){
            refracted = false;
        }
        else{
            refractDir = {(ray.direction[0] * eta) + (N[0] * (eta * cos_i - sqrt(k))),
            (ray.direction[1] * eta) + (N[1] * (eta * cos_i - sqrt(k))),
            (ray.direction[2] * eta) + (N[2] * (eta * cos_i - sqrt(k)))};

            float refractDir_len = sqrt(pow(refractDir[0], 2) + pow(refractDir[1], 2) + pow(refractDir[2], 2));

            refractDir[0] /= refractDir_len;
            refractDir[1] /= refractDir_len;
            refractDir[2] /= refractDir_len;

            refracted = true;
        }

        if (refracted){

            vector<float> refractedOrigin = {hs.hitPoint[0] + (refractDir[0] * EPS),
            hs.hitPoint[1] + (refractDir[1] * EPS),
            hs.hitPoint[2] + (refractDir[2] * EPS)};

            Ray refractRay;
            refractRay.origin = refractedOrigin;
            refractRay.direction = refractDir;


            float rtMin = 0.0f, rtMax = FLT_MAX;
            vector<float> refractColour = {0.0f, 0.0f, 0.0f};
            if (bvh->intersect(refractRay, rtMin, rtMax)) {
                refractColour = reflectRefract(refractRay, bvh, cam, lights, depth + 1);
            } else {
                refractColour = background;
            }

            float cosTheta = fabs(V[0]*N[0] + V[1]*N[1] + V[2]*N[2]);
            float r0 = pow(((n1 - n2) / (n1 + n2)), 2);

            float fresnel = r0 + (1.0f - r0) * pow(1.0f - cosTheta, 5);
            float weight = transparency * (1.0f - fresnel);

            colour[0] += (refractColour[0] * weight);
            colour[1] += (refractColour[1] * weight);
            colour[2] += (refractColour[2] * weight);
        }
    }

    colour[0] = min(1.0f, max(0.0f, colour[0]));
    colour[1] = min(1.0f, max(0.0f, colour[1]));
    colour[2] = min(1.0f, max(0.0f, colour[2]));

    return colour;
}



//----- *Raytracer* -----

int main(){

    vector<Camera> cams = Camera::parseCameraDataFromJson();
    vector<PointLight> lights = PointLight::parseLightDataFromJson();
    vector<Sphere> spheres = Sphere::parseSphereDataFromJson();
    vector<Cube> cubes = Cube::parseCubeDataFromJson();
    vector<Plane> planes = Plane::parsePlaneDataFromJson();

    int totalObjects = spheres.size() + cubes.size() + planes.size();

    bool intersect;
    Ray ray;
    vector<Ray> rays;

    Image img("blank_1920x1080.ppm");
    img.readImage();

    cout << "---Testing for " << totalObjects <<" object(s)---\n";


    /*Accelerated ray caster*/
    //Build BVH
    BVHNode* bvh = BVHNode::buildBVH(planes, cubes, spheres);

    uint64_t start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    //For every pixel
    for (int y = 0; y < cams[0].res_y; ++y) {
        for (int x = 0; x < cams[0].res_x; ++x) {
            intersect = false;

            ray = cams[0].convertPixelToRay(x, y);
            rays.push_back(ray);

            //Traverse BVH
            float tMin = 0.0f;
            float tMax = FLT_MAX;
            if (bvh->intersect(ray, tMin, tMax)) {
                intersect = true;
                HitStructure hs = ray.hs.back();

                //vector<float> colour = blinnPhong(hs, cams[0], lights);

                vector<float> colour = reflectRefract(ray, bvh, cams[0], lights);

                int r = 255 * colour[0];
                int g = 255 * colour[1];
                int b = 255 * colour[2];

                //cout << x << ", " << y << "\n";

                img.modifyPixel(x, y,  r, g, b);
            }
        }
    }
    uint64_t end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();

    img.writeImage("output.ppm");

    cout << "Time taken for BVH RT (accelerated): " << (end - start) / 1000.0f << "s\n";

    delete bvh;

    /*
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
    */

    return 0;
}



