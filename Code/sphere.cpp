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

using namespace std;

//Overarching function for parsing sphere data from JSON
vector<Sphere> Sphere::parseSphereDataFromJson() {
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
        string locationStr = getJSONObject(sphereDataStr, "\"location\"");
        newSphere.location[0] = getFloat(locationStr, "\"x\"");
        newSphere.location[1] = getFloat(locationStr, "\"y\"");
        newSphere.location[2] = getFloat(locationStr, "\"z\"");

        //Parse radius
        newSphere.radius = getFloat(sphereDataStr, "\"radius\"");

        //Parse material data
        string materialStr = getJSONObject(sphereDataStr, "\"material\"");
        newSphere.diffuse[0] = getFloat(materialStr, "\"r\"");
        newSphere.diffuse[1] = getFloat(materialStr, "\"g\"");
        newSphere.diffuse[2] = getFloat(materialStr, "\"b\"");

        string specStr = getJSONObject(sphereDataStr, "\"specular\"");
        newSphere.specular[0] = getFloat(materialStr, "\"r\"");
        newSphere.specular[1] = getFloat(materialStr, "\"g\"");
        newSphere.specular[2] = getFloat(materialStr, "\"b\"");

        newSphere.shininess = getFloat(sphereDataStr, "\"shininess\"");

        newSphere.transparency = getFloat(sphereDataStr, "\"transparency\"");

        newSphere.ior = getFloat(sphereDataStr, "\"ior\"");

        spheres.push_back(newSphere);
    }

    return spheres;
}


bool Sphere::intersect(const Ray& ray, HitStructure& hs){
    
    //Vector from ray origin to sphere centre
    float l[3] = {location[0] - ray.origin[0], location[1] - ray.origin[1], location[2] - ray.origin[2]};
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

    hs.hitPoint = intersectPoint;
    hs.normal = {(intersectPoint[0] - location[0]) / radius, (intersectPoint[1] - location[1]) / radius, (intersectPoint[2] - location[2]) / radius};
    hs.rayDistance = sqrt(pow(ray.origin[0] - intersectPoint[0], 2) + pow(ray.origin[1] - intersectPoint[1], 2) + pow(ray.origin[2] - intersectPoint[2], 2 ));
    hs.diffuse = {diffuse[0], diffuse[1], diffuse[2]};
    hs.specular = {specular[0], specular[1], specular[2]};
    hs.shininess = shininess;
    hs.transparency = transparency;
    hs.ior = ior;

    return true;
}

//Get AABB of the sphere
AABB Sphere::getAABB() const {
    vector<float> min = {location[0] - radius, location[1] - radius, location[2] - radius};
    vector<float> max = {location[0] + radius, location[1] + radius, location[2] + radius};
    return AABB(min, max);
}