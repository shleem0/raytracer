#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "shape.h"
#include "sphere.h"
#include "ray.h"
#include "hit_struct.h"
#include "raytracer.h"
#include "image.h"
#include "config.h"

using namespace std;

//Overarching function for parsing sphere data from JSON
vector<Sphere> Sphere::parseSphereDataFromJson(Config config) {
    //Open JSON file
    ifstream jsonFile("..\\ASCII\\scene.json");
    if (!jsonFile) {
        cerr << "Error opening JSON file." << endl;
        return {};
    }

    //Read file into string
    stringstream buffer;
    buffer << jsonFile.rdbuf();
    string jsonStr = buffer.str();

    //Get spheres array as string
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string spheresArrayStr = getJSONArray(propertiesStr, "\"spheres\"");

    //Remove outer brackets
    if (spheresArrayStr.front() == '[' && spheresArrayStr.back() == ']') {
        spheresArrayStr = spheresArrayStr.substr(1, spheresArrayStr.size() - 2);
    }

    //Split array into sphere objects
    vector<string> sphereObjects;
    int braces = 0;
    size_t start = 0;
    for (size_t i = 0; i < spheresArrayStr.size(); ++i) {
        if (spheresArrayStr[i] == '{') {
            if (braces == 0) start = i;
            braces++;
        } else if (spheresArrayStr[i] == '}') {
            braces--;
            if (braces == 0) {
                sphereObjects.push_back(spheresArrayStr.substr(start, i - start + 1));
            }
        }
    }

    if (sphereObjects.empty()){
        cerr << "No spheres." << endl;
        return {};
    }

    vector<Sphere> spheres;

    for (int i = 0; i < sphereObjects.size(); i++){

        Sphere newSphere;

        //Extract the current sphere
        string sphereDataStr = sphereObjects[i];

        //Parse location
        string locationStr = getJSONObject(sphereDataStr, "\"start_location\"");
        newSphere.startLocation = {getFloat(locationStr, "\"x\""),
        getFloat(locationStr, "\"y\""),
        getFloat(locationStr, "\"z\"")};
        
        locationStr = getJSONObject(sphereDataStr, "\"end_location\"");
        newSphere.endLocation = {getFloat(locationStr, "\"x\""),
        getFloat(locationStr, "\"y\""),
        getFloat(locationStr, "\"z\"")};

        //Parse radius
        newSphere.radius = getFloat(sphereDataStr, "\"radius\"");

        //Parse material data
        string materialStr = getJSONObject(sphereDataStr, "\"material\"");
        
        string diffStr = getJSONObject(materialStr, "\"diffuse\"");
        newSphere.diffuse = {getFloat(diffStr, "\"r\""),
        getFloat(diffStr, "\"g\""),
        getFloat(diffStr, "\"b\"")};

        string specStr = getJSONObject(materialStr, "\"specular\"");
        newSphere.specular = {getFloat(specStr, "\"r\""),
        getFloat(specStr, "\"g\""),
        getFloat(specStr, "\"b\"")};

        newSphere.shininess = getFloat(materialStr, "\"shininess\"");

        newSphere.transparency = getFloat(materialStr, "\"transparency\"");

        newSphere.ior = getFloat(materialStr, "\"ior\"");

        //Parse texture file
        newSphere.texture = getString(materialStr, "\"texture\"");
        if (newSphere.texture != "" && config.textures){
            newSphere.hasTex = true;
            cout << "Sphere " << i << " texture: " << newSphere.texture << "\n";
        }

        spheres.push_back(newSphere);
    }

    return spheres;
}


bool Sphere::intersect(const Ray& ray, HitStructure& hs, Config config){

    vector<float> location;
    if (config.motionBlur){
        location = this->positionAt(ray.time);
    }
    else{
        location = startLocation;
    }
    
    //Vector from ray origin to sphere centre
    vector<float> l = {location[0] - ray.origin[0], location[1] - ray.origin[1], location[2] - ray.origin[2]};
    //Get t closest approach - dot product gives amount of projection of L onto ray's direction vector
    float tca = l[0] * ray.direction[0] + l[1] * ray.direction[1] + l[2] * ray.direction[2];
    
    if (tca < 0){ //If sphere's centre is behind ray origin
        return false; //No intersect
    }

    float d2 = l[0]*l[0] + l[1]*l[1] + l[2]*l[2] - tca*tca;

    if (d2 > radius*radius){ //If ray does not pass through radius
        return false; //No intersectiom
    }

    float thc = sqrt(radius*radius - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;

    float t = (t0 > 0) ? t0 : t1;
    if (t < 0) return false;

    vector<float> intersectPoint = {ray.origin[0] + t * ray.direction[0], ray.origin[1] + t * ray.direction[1], ray.origin[2] + t * ray.direction[2]};

    if (hasTex){
        //Get hit point on texture for texture mapping
        vector<float> normal = Raytracer::normalise(Raytracer::sub_vec(intersectPoint, location));
        float u = 0.5f + atan2(normal[2], normal[0]) / (2 * M_PI);
        float v = 0.5f - asin(normal[1]) / M_PI;

        hs.u = u;
        hs.v = v;
        hs.textureFile = texture;
        hs.hasTex = true;
    }

    hs.hitPoint = intersectPoint;
    hs.normal = {(intersectPoint[0] - location[0]) / radius, (intersectPoint[1] - location[1]) / radius, (intersectPoint[2] - location[2]) / radius};
    hs.rayDistance = sqrt(pow(ray.origin[0] - intersectPoint[0], 2) + pow(ray.origin[1] - intersectPoint[1], 2) + pow(ray.origin[2] - intersectPoint[2], 2 ));
    hs.diffuse = {diffuse[0], diffuse[1], diffuse[2]};
    hs.specular = {specular[0], specular[1], specular[2]};
    hs.shininess = shininess;
    hs.transparency = transparency;
    hs.ior = ior;
    if (ray.time){
        hs.time = ray.time;
    }
    return true;
}

//Get AABB of the sphere
AABB Sphere::getAABB() const {
    vector<float> startMin = {startLocation[0] - radius, startLocation[1] - radius, startLocation[2] - radius};
    vector<float> startMax = {startLocation[0] + radius, startLocation[1] + radius, startLocation[2] + radius};
    AABB start(startMin, startMax);

    vector<float> endMin = {endLocation[0] - radius, endLocation[1] - radius, endLocation[2] - radius};
    vector<float> endMax = {endLocation[0] + radius, endLocation[1] + radius, endLocation[2] + radius};
    AABB end(endMin, endMax);

    AABB fullAABB = AABB::unionOf(start, end);
    return fullAABB;
}