#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include "shape.h"
#include "cube.h"
#include "ray.h"
#include "hit_struct.h"
#include "aabb.h"
#include "image.h"

using namespace std;

//Overarching function for parsing cube data from JSON
vector<Cube> Cube::parseCubeDataFromJson() {
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

    //Get cubes array as string
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string cubesArrayStr = getJSONArray(propertiesStr, "\"cubes\"");

    //Remove outer brackets
    if (cubesArrayStr.front() == '[' && cubesArrayStr.back() == ']') {
        cubesArrayStr = cubesArrayStr.substr(1, cubesArrayStr.size() - 2);
    }

    //Split array into cube objects
    vector<string> cubeObjects;
    int braces = 0;
    size_t start = 0;
    for (size_t i = 0; i < cubesArrayStr.size(); ++i) {
        if (cubesArrayStr[i] == '{') {
            if (braces == 0) start = i;
            braces++;
        } else if (cubesArrayStr[i] == '}') {
            braces--;
            if (braces == 0) {
                cubeObjects.push_back(cubesArrayStr.substr(start, i - start + 1));
            }
        }
    }

    if (cubeObjects.empty()){
        cerr << "No cubes." << endl;
        return {};
    }

    vector<Cube> cubes;

    for (int i = 0; i < cubeObjects.size(); i++){

        Cube newCube;

        //Extract the current cube
        string cubeDataStr = cubeObjects[i];

        //Parse location
        string locationStr = getJSONObject(cubeDataStr, "\"translation\"");
        newCube.location = {getFloat(locationStr, "\"x\""),
        getFloat(locationStr, "\"y\""),
        getFloat(locationStr, "\"z\"")};

        //Parse rotation
        string rotationStr = getJSONObject(cubeDataStr, "\"rotation\"");
        newCube.rotation = {getFloat(rotationStr, "\"x\""),
        getFloat(rotationStr, "\"y\""),
        getFloat(rotationStr, "\"z\"")};

        //Parse 1D scale
        newCube.scale = getFloat(cubeDataStr, "\"scale\"");


        //Parse material data
        string materialStr = getJSONObject(cubeDataStr, "\"material\"");

        string diffStr = getJSONObject(materialStr, "\"diffuse\"");
        newCube.diffuse = {getFloat(diffStr, "\"r\""),
        getFloat(diffStr, "\"g\""),
        getFloat(diffStr, "\"b\"")};

        string specStr = getJSONObject(materialStr, "\"specular\"");
        newCube.specular = {getFloat(specStr, "\"r\""),
        getFloat(specStr, "\"g\""),
        getFloat(specStr, "\"b\"")};

        newCube.shininess = getFloat(materialStr, "\"shininess\"");

        newCube.transparency = getFloat(materialStr, "\"transparency\"");

        newCube.ior = getFloat(materialStr, "\"ior\"");

        newCube.texture = getString(materialStr, "\"texture\"");
        if (newCube.texture != ""){
            cout << "Cube " << i << " texture: " << newCube.texture << "\n";
            newCube.hasTex = true;
        }

        cubes.push_back(newCube);
    }

    return cubes;
}



bool Cube::intersect(const Ray& ray, HitStructure& hs){
    //Transform ray into local space
    Ray local = ray;

    //Translate
    local.origin[0] -= location[0];
    local.origin[1] -= location[1];
    local.origin[2] -= location[2];

    //Inverse rotate
    vector<float> o = {local.origin[0], local.origin[1], local.origin[2]};
    vector<float> d = {local.direction[0]/scale, local.direction[1]/scale, local.direction[2]/scale};
    rotateXYZInverse(o, rotation);
    rotateXYZInverse(d, rotation);
    local.origin = {o[0]/scale, o[1]/scale, o[2]/scale};

    local.direction = d;

    //Intersect with local AABB (-1 to 1)
    float tmin = -INFINITY, tmax = INFINITY;
    for (int i = 0; i < 3; i++) {
        if (fabs(local.direction[i]) < 1e-6f) {
            if (local.origin[i] < -1 || local.origin[i] > 1) return false;
        } else {
            float t1 = (-1 - local.origin[i]) / local.direction[i];
            float t2 = ( 1 - local.origin[i]) / local.direction[i];
            if (t1 > t2) swap(t1, t2);
            tmin = max(tmin, t1);
            tmax = min(tmax, t2);
            if (tmin > tmax) return false;
        }
    }

    float t = (tmin > 0) ? tmin : tmax;
    if (t < 0) return false;

    //Compute hit point and normal (local)
    vector<float> hitLocal= {
        local.origin[0] + t * local.direction[0],
        local.origin[1] + t * local.direction[1],
        local.origin[2] + t * local.direction[2]
    };

    vector<float> normalLocal = {0.0f, 0.0f, 0.0f};
    float eps = 1e-5f;
    float absX = fabs(hitLocal[0]);
    float absY = fabs(hitLocal[1]);
    float absZ = fabs(hitLocal[2]);

    if (absX >= absY - eps && absX >= absZ - eps)
        normalLocal[0] = (hitLocal[0] > 0) ? 1.0f : -1.0f;
    else if (absY >= absX - eps && absY >= absZ - eps)
        normalLocal[1] = (hitLocal[1] > 0) ? 1.0f : -1.0f;
    else
        normalLocal[2] = (hitLocal[2] > 0) ? 1.0f : -1.0f;

        
    if (hasTex){
        float u = 0.0f, v = 0.0f;

        // Determine dominant face axis robustly (choose the largest absolute component)
        float ax = fabs(hitLocal[0]);
        float ay = fabs(hitLocal[1]);
        float az = fabs(hitLocal[2]);

        // choose the axis with the maximum absolute value
        float maxA = max(ax, max(ay, az));

        if (ax >= maxA - eps) {
            // X faces
            if (hitLocal[0] > 0.0f) {
                // +X
                u = ( hitLocal[2] + 1.0f ) * 0.5f;
                v = ( hitLocal[1] + 1.0f ) * 0.5f;
            } else {
                // -X
                // flip u so texture orientation is consistent
                u = ( 1.0f - hitLocal[2] ) * 0.5f;
                v = ( hitLocal[1] + 1.0f ) * 0.5f;
            }
        }
        else if (ay >= maxA - eps) {
            // Y faces
            if (hitLocal[1] > 0.0f) {
                // +Y
                u = ( hitLocal[0] + 1.0f ) * 0.5f;
                v = ( hitLocal[2] + 1.0f ) * 0.5f;
            } else {
                // -Y
                u = ( hitLocal[0] + 1.0f ) * 0.5f;
                v = ( 1.0f - hitLocal[2] ) * 0.5f;
            }
        }
        else {
            // Z faces
            if (hitLocal[2] > 0.0f) {
                // +Z
                u = ( hitLocal[0] + 1.0f ) * 0.5f;
                v = ( hitLocal[1] + 1.0f ) * 0.5f;
            } else {
                // -Z
                u = ( 1.0f - hitLocal[0] ) * 0.5f;
                v = ( hitLocal[1] + 1.0f ) * 0.5f;
            }
        }

        // Wrap/normalize UVs into [0,1) safely
        u = fmodf(u, 1.0f);
        v = fmodf(v, 1.0f);
        if (u < 0.0f) u += 1.0f;
        if (v < 0.0f) v += 1.0f;

        // (optional) clamp tiny numerical issues
        if (u < 0.0f) u = 0.0f;
        if (u > 1.0f) u = 1.0f;
        if (v < 0.0f) v = 0.0f;
        if (v > 1.0f) v = 1.0f;

        hs.u = u;
        hs.v = v;

        hs.hasTex = true;
        hs.textureFile = texture;
    }



    //Transform back to world space
    vector<float> hitWorld = {
        hitLocal[0] * scale,
        hitLocal[1] * scale,
        hitLocal[2] * scale
    };
    rotateXYZ(hitWorld, rotation);
    hitWorld[0] += location[0];
    hitWorld[1] += location[1];
    hitWorld[2] += location[2];

    rotateXYZ(normalLocal, rotation);
    float len = sqrt(normalLocal[0]*normalLocal[0] +
                          normalLocal[1]*normalLocal[1] +
                          normalLocal[2]*normalLocal[2]);
    for (int i = 0; i < 3; i++){
        normalLocal[i] /= len;}

    //Write result
    hs.hitPoint = {hitWorld[0], hitWorld[1], hitWorld[2]};
    hs.normal = normalLocal;
    hs.rayDistance = t;
    hs.diffuse = {diffuse[0], diffuse[1], diffuse[2]};
    hs.specular = {specular[0], specular[1], specular[2]};
    hs.shininess = shininess;
    hs.transparency = transparency;
    hs.ior = ior;

    return true;
}


//Helper functions

void Cube::rotateXYZ(vector<float>& v, vector<float>& rot) {
    float x = v[0], y = v[1], z = v[2];

    //X rotation
    float cx = cos(rot[0]), sx = sin(rot[0]);
    float y1 = y * cx - z * sx;
    float z1 = y * sx + z * cx;
    y = y1; z = z1;

    //Y rotation
    float cy = cos(rot[1]), sy = sin(rot[1]);
    float x2 = x * cy + z * sy;
    float z2 = -x * sy + z * cy;
    x = x2; z = z2;

    //Z rotation
    float cz = cos(rot[2]), sz = sin(rot[2]);
    float x3 = x * cz - y * sz;
    float y3 = x * sz + y * cz;

    v[0] = x3; v[1] = y3; v[2] = z2;
}

void Cube::rotateXYZInverse(vector<float>& v, vector<float>& rot) {
    vector<float> inv = {-rot[0], -rot[1], -rot[2]};
    rotateXYZ(v, inv);
}

AABB Cube::getAABB() const {
    vector<float> min = {location[0] - scale, location[1] - scale, location[2] - scale};
    vector<float> max = {location[0] + scale, location[1] + scale, location[2] + scale};
    return AABB(min, max);
}
