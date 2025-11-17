#include "raytracer.h"
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
#include <algorithm>
#include <utility>
using namespace std;
using namespace std::chrono;



//----- *Raytracer* -----

vector<pair<string, Image>> textureCache;

int main(){

    //Parse object data from JSON
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


    cout << "---Raytracing for " << totalObjects <<" object(s)---\n";

    /*Accelerated ray caster*/
    //Build BVH
    BVHNode* bvh = BVHNode::buildBVH(planes, cubes, spheres);

    int samplesPerPixel = 4;
    float invSamples = 1.0f / samplesPerPixel;

    uint64_t start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    //For every pixel
    for (int y = 0; y < cams[0].res_y; ++y) {
        for (int x = 0; x < cams[0].res_x; ++x) {

            intersect = false;
            vector<float> pixelColour = {0.0f, 0.0f, 0.0f};


            //Random (Monte Carlo) samples for anti-aliasing
            for (int s = 0; s < samplesPerPixel; s++){

                //Get random position variation within pixel
                float u = x + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                float v = y + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

                //Get the related ray
                ray = cams[0].convertPixelToRay(u, v);
                rays.push_back(ray);

                //Traverse BVH
                float tMin = 0.0f;
                float tMax = FLT_MAX;
                ray.hs.clear();
                if (bvh->intersect(ray, tMin, tMax)) {
                    intersect = true;
                    HitStructure hs = ray.hs.back();

                    //Get Blinn-Phong, reflection and refraction components of pixel lighting
                    vector<float> colour = Raytracer::reflectRefract(ray, bvh, cams[0], lights);
                    pixelColour = Raytracer::add_vec(pixelColour, colour);
                }
            }

            //Average the colour for AA
            pixelColour = Raytracer::mul_vec(pixelColour, invSamples);

            int r = 255 * pixelColour[0];
            int g = 255 * pixelColour[1];
            int b = 255 * pixelColour[2];

            //Edit the template image
            img.modifyPixel(x, y,  r, g, b);
        }
    }
    uint64_t end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    cout << "Time taken for BVH RT (accelerated): " << (end - start) / 1000.0f << "s\n";

    img.writeImage("output.ppm");

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



//Blinn-Phong shading for a given ray
vector<float> Raytracer::blinnPhong(const HitStructure& hs, const Camera& cam, const vector<PointLight> lights){

    float ka = 0.5f; //Ambient reflect
    float kd = 0.2126f * hs.diffuse[0] + 0.7152f * hs.diffuse[1] + 0.0722f * hs.diffuse[2]; //Diffuse reflect
    float ks = 0.2126f * hs.specular[0] + 0.7152f * hs.specular[1] + 0.0722f * hs.specular[2]; //Specular reflect
    float shiny = hs.shininess;
    vector<float> matColour;

    //Gather colour from texture
    if (hs.hasTex){

        Image* texPtr = nullptr;
        //Check if in texture cache
        for (auto& pair : textureCache){
            if (pair.first == hs.textureFile){
                texPtr = &pair.second;
                break;
            }
        }
        if (!texPtr){
            textureCache.push_back({hs.textureFile, Image(hs.textureFile)});
            texPtr = &textureCache.back().second;
            texPtr->readImage();
        }

        float u = fmod(hs.u, 1.0f); if (hs.u < 0) u += 1.0f;
        float v = fmod(hs.v, 1.0f); if (hs.v < 0) v += 1.0f;

        int x = std::min(int(u * (texPtr->width - 1)), texPtr->width - 1);
        int y = std::min(int((1 - v) * (texPtr->height - 1)), texPtr->height - 1);

        //Extract from texture
        matColour = {texPtr->pixelValues[y][x][0] / 255.0f,
        texPtr->pixelValues[y][x][1] / 255.0f,
        texPtr->pixelValues[y][x][2] / 255.0f};
    }
    else{
        matColour = hs.diffuse;
    }
    vector<float> N = hs.normal;

    //View vector
    vector<float> V = {cam.location[0] - hs.hitPoint[0], cam.location[1] - hs.hitPoint[1], cam.location[2] - hs.hitPoint[2]};
    V = normalise(V);


    //Ambient colour
    vector<float> colour = {ka * matColour[0], ka * matColour[1], ka * matColour[2]};

    for (const auto& light : lights){

        //Get vector to each light from the hit point
        vector<float> L = {light.location[0] - hs.hitPoint[0], light.location[1] - hs.hitPoint[1], light.location[2] - hs.hitPoint[2]};
        float L_len = sqrt(pow(L[0], 2) + pow(L[1], 2) + pow(L[2], 2));
        L = normalise(L);

        float irr = light.rad_intensity / (1.0f + pow(L_len, 2));

        vector<float> H = add_vec(L, V);
        H = normalise(H);

        float N_L = max(0.0f, dotProd(N, L));
        float N_H = max(0.0f, dotProd(N, H));

        //Adjust colour using the irradiance of the light and specularity of the material
        colour[0] += matColour[0] * N_L * irr + 0.5f * hs.specular[0] * pow(N_H, shiny) * irr;
        colour[1] += matColour[1] * N_L * irr + 0.5f * hs.specular[1] * pow(N_H, shiny) * irr;
        colour[2] += matColour[2] * N_L * irr + 0.5f * hs.specular[2] * pow(N_H, shiny) * irr;
    }

    colour[0] = min(1.0f, max(colour[0], 0.0f));
    colour[1] = min(1.0f, max(colour[1], 0.0f));
    colour[2] = min(1.0f, max(colour[2], 0.0f));

    return colour;
}



//Recursive calculation of reflection and refraction rays at a given hitpoint
vector<float> Raytracer::reflectRefract(Ray ray, BVHNode* bvh, const Camera& cam, const vector<PointLight>& lights, int depth){

    const int MAX_DEPTH = 5;
    const float EPS = 1e-4f;
    const vector<float> background = {0.0f, 0.0f, 0.0f};
    float tMin = 0.0f, tMax = FLT_MAX;

    //Stop recursion at max depth, or if no intersection occurs
    if (depth > MAX_DEPTH || !bvh->intersect(ray, tMin, tMax)){
        return background;
    }

    HitStructure hs = ray.hs.back();

    vector<float> colour = blinnPhong(hs, cam, lights);

    float reflectivity = ((hs.specular[0] + hs.specular[1] + hs.specular[2]) / 3.0f) * 0.5;
    reflectivity = min(max(reflectivity, 0.0f), 1.0f);

    float ior = hs.ior > 0.0f ? hs.ior : 1.0f;
    float transparency = hs.transparency;

    vector<float> V = sub_vec(cam.location, hs.hitPoint);
    V = normalise(V);

    vector<float> N = hs.normal;

    float N_rayDir = dotProd(ray.direction, N);
    if (N_rayDir > 0.0f) {
        //Flip normal so N faces opposite to ray.direction
        N = {-N[0], -N[1], -N[2]};
    }


    //Reflection
    vector<float> reflectColour = {0.0f, 0.0f, 0.0f};
    if (reflectivity > 0){

        float incoming_N = dotProd(ray.direction, N);
    
        vector<float> reflectDir = {ray.direction[0] - 2.0f * incoming_N * N[0],
        ray.direction[1] - 2.0f * incoming_N * N[1],
        ray.direction[2] - 2.0f * incoming_N * N[2]
        };

        reflectDir = normalise(reflectDir);

        vector<float> reflectEPSDir = mul_vec(N, EPS);
        vector<float> reflectOrigin = add_vec(hs.hitPoint, reflectEPSDir);

        Ray reflectRay;
        reflectRay.origin = reflectOrigin;
        reflectRay.direction = reflectDir;
        reflectRay.hs.clear();

        float rtMin = 0.0f, rtMax = FLT_MAX;
        ray.hs.clear();
        if (bvh->intersect(reflectRay, rtMin, rtMax)) {
            reflectColour = reflectRefract(reflectRay, bvh, cam, lights, depth + 1);
        } 
        else {
            reflectColour = background;
        }
    }

    reflectColour = mul_vec(reflectColour, reflectivity);
    colour = add_vec(colour, reflectColour);
    

    //Refraction
    bool refracted = false;
    if (transparency > 0){
        float n1 = 1.0f;
        float n2 = ior;
        float cos_i = dotProd(ray.direction, N);

        if (cos_i > 0){
            swap(n1, n2);
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

            refractDir = normalise(refractDir);

            refracted = true;
        }

        if (refracted){

            vector<float> refractEPSDir = mul_vec(refractDir, EPS);
            vector<float> refractedOrigin = add_vec(hs.hitPoint, refractEPSDir);

            Ray refractRay;
            refractRay.origin = refractedOrigin;
            refractRay.direction = refractDir;


            float rtMin = 0.0f, rtMax = FLT_MAX;
            vector<float> refractColour = {0.0f, 0.0f, 0.0f};
            ray.hs.clear();
            if (bvh->intersect(refractRay, rtMin, rtMax)) {
                refractColour = reflectRefract(refractRay, bvh, cam, lights, depth + 1);
            } else {
                refractColour = background;
            }

            float cosTheta = fabs(dotProd(V, N));
            float r0 = pow(((n1 - n2) / (n1 + n2)), 2);

            float fresnel = r0 + (1.0f - r0) * pow(1.0f - cosTheta, 5);
            float weight = transparency * (1.0f - fresnel);


            refractColour = mul_vec(refractColour, weight);
            colour = add_vec(colour, refractColour);
        }
    }

    colour[0] = min(1.0f, max(0.0f, colour[0]));
    colour[1] = min(1.0f, max(0.0f, colour[1]));
    colour[2] = min(1.0f, max(0.0f, colour[2]));

    return colour;
}



//----- Helper functions -----

//Normalise vector
vector<float> Raytracer::normalise(vector<float> vec){
    float vec_len = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));

    vec[0] /= vec_len;
    vec[1] /= vec_len;
    vec[2] /= vec_len;

    return vec;
}

//Get dot product of 2 vectors
float Raytracer::dotProd(vector<float> v1, vector<float> v2){
    float result = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

    return result;
}

//Add 2 vectors
vector<float> Raytracer::add_vec(vector<float> v1, vector<float> v2){
    vector<float> result = {v1[0] + v2[0],
    v1[1] + v2[1],
    v1[2] + v2[2]};

    return result;
}


//Subtract 2 vectors
vector<float> Raytracer::sub_vec(vector<float> v1, vector<float> v2){
    vector<float> result = {v1[0] - v2[0],
    v1[1] - v2[1],
    v1[2] - v2[2]};

    return result;
}


//Multiply a vector by a scalar
vector<float> Raytracer::mul_vec(vector<float> vec, float mul){
    vector<float> result = {vec[0] * mul,
    vec[1] * mul,
    vec[2] * mul};

    return result;
}

