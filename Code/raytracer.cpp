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
#include "config.h"
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



//----- *Main Raytracer* -----
vector<pair<string, Image>> textureCache;

int main(int argc, char* argv[]){

    Config config = Config::parseArgs(argc, argv);

    cout << "---Config---"
        << "\nBVH: " << config.bvh
        << "\nAnti-aliasing: " << config.antiAliasing
        << "\n  AA Samples: " << config.AAsamples
        << "\nReflections: " << config.reflectRefractDepth
        << "\nTexture Mapping: " << config.textures
        << "\n-Distributed effects-"
        << "\n  Glossy reflections: " << config.glossyReflect
        << "\n  GR samples: " << config.GRsamples
        << "\n  Soft shadows: " << config.softShadows
        << "\n  SS samples: " << config.SSsamples
        << "\n-Lens effects-"
        << "\n  Motion blur: " << config.motionBlur
        << "\n  MB samples: " << config.mbSamples
        << "\n  Depth of Field: " << config.dof
        << "\n  DoF samples: " << config.dofSamples << "\n\n";




    //Parse object data from JSON
    vector<Camera> cams = Camera::parseCameraDataFromJson(config);
    vector<PointLight> lights = PointLight::parseLightDataFromJson();
    vector<Sphere> spheres = Sphere::parseSphereDataFromJson(config);
    vector<Cube> cubes = Cube::parseCubeDataFromJson(config);
    vector<Plane> planes = Plane::parsePlaneDataFromJson(config);

    int totalObjects = spheres.size() + cubes.size() + planes.size();

    Image img("blank_1920x1080.ppm");
    img.readImage();

    cout << "\n---Raytracing for " << totalObjects <<" object(s)---\n";

    BVHNode* bvh;
    Ray ray;

    //Build BVH
    if (config.bvh){
        bvh = BVHNode::buildBVH(planes, cubes, spheres, config);
        cout << "Running BVH accelerated RT...\n";
    }
    else{
        bvh = nullptr;
        cout << "Running unaccelerated RT...\n";
    }

    int AAsamplesPerPixel;
    if (config.antiAliasing){
        AAsamplesPerPixel = config.AAsamples;
    }
    else{
        AAsamplesPerPixel = 1;
    }

    int dofSamplesPerPixel;
    if (config.dof){
        dofSamplesPerPixel = config.dofSamples;
    }
    else{
        dofSamplesPerPixel = 1;
    }

    int mbSamplesPerPixel;
    if (config.motionBlur){
        mbSamplesPerPixel = config.mbSamples;
    }
    else{
        mbSamplesPerPixel = 1;
    }

    int percentDone = 0;
    int pixelsDone = 0;
    int totalPixels = static_cast<int>(cams[0].res_y * cams[0].res_x);

    uint64_t start = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    //For every pixel
    for (int y = 0; y < cams[0].res_y; ++y) {
        for (int x = 0; x < cams[0].res_x; ++x) {

            vector<float> pixelColour = {0.0f, 0.0f, 0.0f};

            float u;
            float v;

            //Random (Monte Carlo) samples for anti-aliasing
            for (int s = 0; s < AAsamplesPerPixel; s++){


                if (config.antiAliasing){
                //Get random position variation within pixel
                    u = x + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                    v = y + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                }
                else{
                    u = x;
                    v = y;
                }

                for (int t = 0; t < dofSamplesPerPixel; t++){

                    //Get the related ray
                    ray = cams[0].convertPixelToRay(u, v);

                    for (int m = 0; m < mbSamplesPerPixel; m++){

                        ray.hs.clear();
                        ray.time = (m + static_cast<float>(rand())/RAND_MAX) / mbSamplesPerPixel;

                        //Accelerated RT with BVH
                        if (config.bvh){

                            float tMin = 0.0f;
                            float tMax = FLT_MAX;
                            if (bvh->intersect(ray, tMin, tMax, config)) {

                                HitStructure hs = ray.hs.back();

                                //Get Blinn-Phong, reflection and refraction components of pixel lighting
                                vector<float> colour = Raytracer::reflectRefract(config, bvh, ray, cams[0], lights, planes, cubes, spheres, 0, config.reflectRefractDepth);
                                pixelColour = Raytracer::add_vec(pixelColour, colour);
                            }
                        }
                        
                        //Unaccelerated RT
                        else{
                            if (Raytracer::unacceleratedIntersection(ray, planes, spheres, cubes, config)){
                                vector<float> colour = Raytracer::reflectRefract(config, bvh, ray, cams[0], lights, planes, cubes, spheres, 0, config.reflectRefractDepth);
                                pixelColour = Raytracer::add_vec(pixelColour, colour);
                            }
                        }
                    }
                }
            }
    

            //Average the colour for AA
            pixelColour = Raytracer::mul_vec(pixelColour, 1.0f / (AAsamplesPerPixel*dofSamplesPerPixel*mbSamplesPerPixel));

            int r = 255 * pixelColour[0];
            int g = 255 * pixelColour[1];
            int b = 255 * pixelColour[2];

            //Edit the template image
            img.modifyPixel(x, y, r, g, b);
            pixelsDone += 1;

            int newPercent = (pixelsDone * 100) / totalPixels;
            if (newPercent >= percentDone + 10 && newPercent != 100){
                percentDone += 10;
                cout << percentDone << "% done\n";

            }

        }
    }

    //End timer
    uint64_t end = duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    cout << "Finished! Time taken: " << (end - start) / 1000.0f << "s\n";

    img.writeImage(config.outputFile);
    if (bvh) delete bvh;

    return 0;
}

//Unaccelerated object-by-object intersection testing
bool Raytracer::unacceleratedIntersection(Ray& ray, vector<Plane>& planes, vector<Sphere>& spheres, vector<Cube>& cubes, Config config){

    HitStructure hs;
    bool hit = false;
    float closestDist = FLT_MAX;
    bool intersect;

    //Check intersection with planes
    for (int p = 0; p < planes.size(); p++){

        HitStructure temp;
        intersect = planes[p].intersect(ray, temp, config);

        if (intersect && temp.rayDistance < closestDist){

            closestDist = temp.rayDistance;
            hs = temp;
            hit = true;
        }
    }

    //Check intersection with cubes
    for (int c = 0; c < cubes.size(); c++){

        HitStructure temp;
        intersect = cubes[c].intersect(ray, temp, config);

        if (intersect && temp.rayDistance < closestDist){
            closestDist = temp.rayDistance;
            hs = temp;
            hit = true;
        }
    }

    //Check intersection with spheres
    for (int s = 0; s < spheres.size(); s++){

        HitStructure temp;
        intersect = spheres[s].intersect(ray, temp, config);

        if (intersect && temp.rayDistance < closestDist){
            closestDist = temp.rayDistance;
            hs = temp;
            hit = true;
        }
    }

    if (hit) {
        ray.hs.push_back(hs);
        return true;
    }
    else{
        return false;
    }
}




//Blinn-Phong shading at a given point
vector<float> Raytracer::blinnPhong(Config config, BVHNode* bvh, HitStructure& hs, const Camera& cam, vector<PointLight>& lights, vector<Plane>& planes, vector<Cube>& cubes, vector<Sphere>& spheres){

    float ka = 0.25f; //Ambient reflect
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
        }//Otherwise add to cache
        if (!texPtr){
            textureCache.push_back({hs.textureFile, Image(hs.textureFile)});
            texPtr = &textureCache.back().second;
            texPtr->readImage();
        }

        float u = fmod(hs.u, 1.0f); if (hs.u < 0) u += 1.0f;
        float v = fmod(hs.v, 1.0f); if (hs.v < 0) v += 1.0f;

        int x = std::min(int(u * (texPtr->width - 1)), texPtr->width - 1);
        int y = std::min(int((1 - v) * (texPtr->height - 1)), texPtr->height - 1);

        //Extract colour from texture
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

        //Get irradiance at point
        float irr = light.rad_intensity / (1.0f + pow(L_len, 2));

        vector<float> H = Raytracer::add_vec(L, V);
        H = Raytracer::normalise(H);

        float N_L = max(0.0f, dotProd(N, L));
        float N_H = max(0.0f, dotProd(N, H));

        float shadow;

        if (config.softShadows){
            shadow = Raytracer::computeSoftShadows(hs.hitPoint, hs.time, N, light, bvh, planes, cubes, spheres, config);
        }
        else{
            shadow = Raytracer::computeHardShadows(hs.hitPoint, hs.time, N, light, bvh, planes, cubes, spheres, config);
        }
        //Adjust colour using the irradiance of the light and specularity of the material
        colour[0] += shadow * matColour[0] * N_L * irr + 0.5f * hs.specular[0] * pow(N_H, shiny) * irr;
        colour[1] += shadow * matColour[1] * N_L * irr + 0.5f * hs.specular[1] * pow(N_H, shiny) * irr;
        colour[2] += shadow * matColour[2] * N_L * irr + 0.5f * hs.specular[2] * pow(N_H, shiny) * irr;
    }

    colour[0] = min(1.0f, max(colour[0], 0.0f));
    colour[1] = min(1.0f, max(colour[1], 0.0f));
    colour[2] = min(1.0f, max(colour[2], 0.0f));

    return colour;
}



//Recursive calculation of reflection and refraction rays at a given hitpoint
vector<float> Raytracer::reflectRefract(Config config, BVHNode* bvh, Ray& ray, const Camera& cam, vector<PointLight>& lights, vector<Plane>& planes, vector<Cube>& cubes, vector<Sphere>& spheres, int depth, int maxDepth){

    const float EPS = 1e-4f;
    const vector<float> background = {0.0f, 0.0f, 0.0f};
    float tMin = 0.0f, tMax = FLT_MAX;
    bool intersect;

    //If no hits, return background
    if (ray.hs.empty()){
        return background;
    }

    HitStructure hs = ray.hs.back();

    vector<float> colour = blinnPhong(config, bvh, hs, cam, lights, planes, cubes, spheres);

    //Stop reflection/refraction at max depth
    if (depth >= maxDepth){
        return colour;
    }

    //Get reflectivity from average specularity (and clamp)
    float reflectivity = ((hs.specular[0] + hs.specular[1] + hs.specular[2]) / 3.0f) * 0.5;
    reflectivity = min(max(reflectivity, 0.0f), 1.0f);

    float ior = hs.ior > 0.0f ? hs.ior : 1.0f;
    float transparency = hs.transparency;

    //Direction from camera to hitpoint
    vector<float> V = sub_vec(cam.location, hs.hitPoint);
    V = normalise(V);

    vector<float> N = hs.normal;

    float N_rayDir = dotProd(ray.direction, N);
    if (N_rayDir > 0.0f) {
        //Flip normal so N faces opposite to ray.direction
        N = {-N[0], -N[1], -N[2]};
    }

    //Reflection
    vector<float> sampleColour;
    vector<float> reflectColour = {0.0f, 0.0f, 0.0f};
    if (reflectivity > 0){

        float incoming_N = Raytracer::dotProd(ray.direction, N);

        //Do for number of glossy reflection rays
        for (int s = 0; s < config.GRsamples; s++){

            //Ideal reflection vector direction
            vector<float> reflectDir = {ray.direction[0] - 2.0f * incoming_N * N[0],
            ray.direction[1] - 2.0f * incoming_N * N[1],
            ray.direction[2] - 2.0f * incoming_N * N[2]};

            reflectDir = normalise(reflectDir);
            
            if (config.glossyReflect && hs.shininess > 0.0f){

                float glossAngle = max(0.01f, (1.0f - min(1.0f, hs.shininess / 128.0f))) * (M_PI / 6.0f);
                reflectDir = rndConeDirection(reflectDir, glossAngle);
            }

            vector<float> reflectEPSDir = Raytracer::mul_vec(N, EPS);
            vector<float> reflectOrigin = Raytracer::add_vec(hs.hitPoint, reflectEPSDir);

            Ray reflectRay;
            reflectRay.origin = reflectOrigin;
            reflectRay.direction = reflectDir;
            reflectRay.hs.clear();

            float rtMin = 0.0f, rtMax = FLT_MAX;
            if ((bvh && bvh->intersect(reflectRay, rtMin, rtMax, config)) || unacceleratedIntersection(reflectRay, planes, spheres, cubes, config)) {
                sampleColour = reflectRefract(config, bvh, reflectRay, cam, lights, planes, cubes, spheres, depth + 1, maxDepth);
            } 
            else {
                sampleColour = background;
            }

            reflectColour = Raytracer::add_vec(reflectColour, sampleColour);
        }
    }

    //Average over samples for glossy reflection
    reflectColour = Raytracer::mul_vec(reflectColour, 1.0f / config.GRsamples);
    //Scale by reflectivity of material
    reflectColour = Raytracer::mul_vec(reflectColour, reflectivity);
    colour = Raytracer::add_vec(colour, reflectColour);
    

    //Refraction
    bool refracted = false;
    if (transparency > 0){
        float n1 = 1.0f;
        float n2 = ior;
        float cos_i = Raytracer::dotProd(ray.direction, N);

        if (cos_i > 0){
            swap(n1, n2);
        }

        //Check if ray is refracted using refractive index
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

            refractDir = Raytracer::normalise(refractDir);
            refracted = true;
        }

        if (refracted){

            //Get refracted ray
            vector<float> refractEPSDir = Raytracer::mul_vec(refractDir, EPS);
            vector<float> refractedOrigin = Raytracer::add_vec(hs.hitPoint, refractEPSDir);

            Ray refractRay;
            refractRay.origin = refractedOrigin;
            refractRay.direction = refractDir;


            float rtMin = 0.0f, rtMax = FLT_MAX;
            vector<float> refractColour = {0.0f, 0.0f, 0.0f};
            if ((bvh && bvh->intersect(refractRay, rtMin, rtMax, config)) || unacceleratedIntersection(refractRay, planes, spheres, cubes, config)) {
                refractColour = reflectRefract(config, bvh, refractRay, cam, lights, planes, cubes, spheres, depth + 1, maxDepth);
            } else {
                refractColour = background;
            }

            float cosTheta = fabs(Raytracer::dotProd(V, N));
            float r0 = pow(((n1 - n2) / (n1 + n2)), 2);

            float fresnel = r0 + (1.0f - r0) * pow(1.0f - cosTheta, 5);
            float weight = transparency * (1.0f - fresnel);

            refractColour = Raytracer::mul_vec(refractColour, weight);
            colour = Raytracer::add_vec(colour, refractColour);
        }
    }

    colour[0] = min(1.0f, max(0.0f, colour[0]));
    colour[1] = min(1.0f, max(0.0f, colour[1]));
    colour[2] = min(1.0f, max(0.0f, colour[2]));

    return colour;
}

//Binary shadowed/unshadowed calculation for a point
float Raytracer::computeHardShadows(const vector<float>& hitPoint, const float time, const vector<float>& norm, const PointLight& light, BVHNode* bvh, vector<Plane>& planes, vector<Cube>& cubes, vector<Sphere>& spheres, Config config) {
    vector<float> toLight = sub_vec(light.location, hitPoint);
    float lightDist = sqrt(dotProd(toLight, toLight));
    vector<float> L = normalise(toLight);

    Ray shadowRay;
    shadowRay.origin = add_vec(hitPoint, mul_vec(norm, 1e-2f));
    shadowRay.direction = L;
    shadowRay.time = time;

    //If intersection between point and light, there is a shadow
    float tMin = 1e-4f, tMax = lightDist;
    bool hit;
    if (bvh){
        hit = bvh->intersect(shadowRay, tMin, tMax, config);
    }
    else{
        hit = unacceleratedIntersection(shadowRay, planes, spheres, cubes, config);
    }

    return hit ? 0.0f : 1.0f;
}


//Calculate soft shadows by average light blockage over a small area
float Raytracer::computeSoftShadows(const vector<float>& hitPoint, const float time, const vector<float>& norm, const PointLight& light, BVHNode* bvh, vector<Plane>& planes, vector<Cube>& cubes, vector<Sphere>& spheres, Config config){

    int unblocked = 0;

    for (int i = 0; i < config.SSsamples; i++){

        //Slight adjustments over radius to get different light rays
        vector<float> jitter = Raytracer::mul_vec(Raytracer::rndUnitSphere(), config.lightRadius);
        vector<float> lightPos = Raytracer::add_vec(light.location, jitter);

        vector<float> hitToLight = Raytracer::sub_vec(lightPos, hitPoint);

        vector<float> lightDir = Raytracer::normalise(hitToLight);

        Ray shadowRay;
        shadowRay.origin = Raytracer::add_vec(hitPoint, Raytracer::mul_vec(lightDir, 0.001f));
        shadowRay.direction = lightDir;
        shadowRay.time = time;

        float lightDist = sqrt(pow(hitToLight[0], 2) + pow(hitToLight[1], 2) + pow(hitToLight[2], 2));

        bool hit;
        if (bvh){
            float rtMin = 0.0f, rtMax = lightDist;
            hit = bvh->intersect(shadowRay, rtMin, rtMax, config);
        }
        else{
            hit = Raytracer::unacceleratedIntersection(shadowRay, planes, spheres, cubes, config);
        }

        if (!(hit && shadowRay.hs.back().rayDistance < lightDist)){
            unblocked++;
        }
    }

    //Get shadowed ratio
    return static_cast<float>(unblocked) / config.SSsamples;
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

//Get cross product of 2 vectors
vector<float> Raytracer::crossProd(vector<float> v1, vector<float> v2){

    vector<float> result = {v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0]};
    
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

//Sample random point within unit sphere
vector<float> Raytracer::rndUnitSphere(){
    float u = static_cast<float>(rand()) / RAND_MAX;
    float v = static_cast<float>(rand()) / RAND_MAX;
    float theta = u * 2.0f * M_PI;
    float phi = acos(2.0f * v - 1.0f);
    float r = cbrt(static_cast<float>(rand()) / RAND_MAX);

    return {r * sin(phi) * cos(theta),
    r * sin(phi) * sin(theta),
    r * cos(phi)};
}

//Sample random direction from cone
vector<float> Raytracer::rndConeDirection(const vector<float>& dir, float angleRad) {
    vector<float> randVec = rndUnitSphere(); // unit vector in sphere
    float cosTheta = cos(angleRad);
    float z = cosTheta + (1.0f - cosTheta) * ((rand() % 1000) / 1000.0f);
    float phi = 2.0f * M_PI * ((rand() % 1000) / 1000.0f);
    float r = sqrt(1.0f - z*z);
    randVec = {r * cos(phi), r * sin(phi), z};

    // Build an orthonormal basis with dir as Z
    vector<float> w = dir;
    vector<float> u = Raytracer::normalise(Raytracer::crossProd({0.0f,1.0f,0.0f}, w));
    vector<float> v = Raytracer::crossProd(w, u);

    // Convert randVec from local space to world space
    vector<float> sampleDir = {
        randVec[0]*u[0] + randVec[1]*v[0] + randVec[2]*w[0],
        randVec[0]*u[1] + randVec[1]*v[1] + randVec[2]*w[1],
        randVec[0]*u[2] + randVec[1]*v[2] + randVec[2]*w[2]
    };

    return Raytracer::normalise(sampleDir);
}

//Linear interpolation
float Raytracer::lerp(float a, float b, float t) {
    return a + t * (b - a);
}




