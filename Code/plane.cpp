#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include "shape.h"
#include "plane.h"
#include "ray.h"
#include "hit_struct.h"

using namespace std;

vector<Plane> Plane::parsePlaneDataFromJson() {
    // Open JSON file
    ifstream jsonFile("..\\ASCII\\scene.json");
    if (!jsonFile) {
        cerr << "Error opening JSON file." << endl;
        return {};
    }

    // Read entire file into string
    stringstream buffer;
    buffer << jsonFile.rdbuf();
    string jsonStr = buffer.str();

    // Get planes array as string (includes [ ... ])
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string planesArrayStr = getJSONArray(propertiesStr, "\"planes\"");

    // Remove outer brackets
    if (planesArrayStr.front() == '[' && planesArrayStr.back() == ']') {
        planesArrayStr = planesArrayStr.substr(1, planesArrayStr.size() - 2);
    }

    // Split array into full plane objects
    vector<string> planeObjects;
    int braces = 0;
    size_t start = 0;
    for (size_t i = 0; i < planesArrayStr.size(); ++i) {
        if (planesArrayStr[i] == '{') {
            if (braces == 0) start = i;
            braces++;
        } else if (planesArrayStr[i] == '}') {
            braces--;
            if (braces == 0) {
                planeObjects.push_back(planesArrayStr.substr(start, i - start + 1));
            }
        }
    }

    if (planeObjects.empty()){
        cerr << "No planes." << endl;
        return {};
    }

    vector<Plane> planes;

    for (int i = 0; i < planeObjects.size(); i++){

        Plane newPlane;

        // Extract the specific plane
        string planeDataStr = planeObjects[i];

        // Parse corners
        string cornersStr = getJSONArray(planeDataStr, "\"corners\"");
        if (cornersStr.front() == '[' && cornersStr.back() == ']') {
            cornersStr = cornersStr.substr(1, cornersStr.size() - 2);
        }

        vector<string> cornerObjects;
        int cornersBraces = 0;
        size_t cornerStart = 0;
        for (size_t i = 0; i < cornersStr.size(); ++i) {
            if (cornersStr[i] == '{') {
                if (cornersBraces == 0) cornerStart = i;
                cornersBraces++;
            } else if (cornersStr[i] == '}') {
                cornersBraces--;
                if (cornersBraces == 0) {
                    cornerObjects.push_back(cornersStr.substr(cornerStart, i - cornerStart + 1));
                }
            }
        }

        for (int i = 0; i < cornerObjects.size(); i++) {


            // Extract the specific corner
            string cornerDataStr = cornerObjects[i];

            float x = getFloat(cornerDataStr, "\"x\"");
            float y = getFloat(cornerDataStr, "\"y\"");
            float z = getFloat(cornerDataStr, "\"z\"");

            vector<float> pos = {x, y, z};
            newPlane.vertices.push_back(pos);

        }

        newPlane.calculateNormal();
        planes.push_back(newPlane);
    }

    return planes;
}

bool Plane::intersect(const Ray& ray, HitStructure& hs){
    float denominator = normal[0]*ray.direction[0] + normal[1]*ray.direction[1] + normal[2]*ray.direction[2];
    if (abs(denominator) < 1e-10){ 
        cout << "parallel\n";
        return false; // Ray is parallel to plane
    }
    
    float t = (normal[0]*(vertices[0][0] - ray.origin[0]) +
            normal[1]*(vertices[0][1] - ray.origin[1]) +
            normal[2]*(vertices[0][2] - ray.origin[2])) / denominator;
    if (t < 0){
        return false; // Intersection is behind the ray
    }

    vector<float> intersectPoint = {ray.origin[0] + t*ray.direction[0], ray.origin[1] + t*ray.direction[1], ray.origin[2] + t*ray.direction[2]}; 

    sortVerticesWinding();
    bool inPolygon = pointInsidePolygon(intersectPoint);

    if (inPolygon){
        hs.hitPoint = intersectPoint;
        hs.rayDistance = sqrt(pow(ray.origin[0] - intersectPoint[0], 2) + pow(ray.origin[1] - intersectPoint[1], 2) + pow(ray.origin[2] - intersectPoint[2], 2 ));
    }

    return inPolygon;
}


/*---Helper functions---*/
bool Plane::pointInsidePolygon(const vector<float>& point){
    //Choose drop axis on dominant normal component
    int dropAxis;
    float nx = fabs(normal[0]);
    float ny = fabs(normal[1]);
    float nz = fabs(normal[2]);
    if (nx > ny && nx > nz)
        dropAxis = 0; // drop X
    else if (ny > nz)
        dropAxis = 1; // drop Y
    else
        dropAxis = 2; // drop Z

    // Project 3D point into 2D (drop one coord)
    auto project2D = [&](const vector<float>& p) -> pair<float, float> {
        switch (dropAxis) {
            case 0:  return {p[1], p[2]}; // drop X
            case 1:  return {p[0], p[2]}; // drop Y
            default: return {p[0], p[1]}; // drop Z
        }
    };

    pair<float, float> P = project2D(point);

    // 2D ray-casting algorithm
    bool inside = false;
    int n = static_cast<int>(vertices.size());
    for (int i = 0, j = n - 1; i < n; j = i++) {
        pair<float, float> Pi = project2D(vertices[i]);
        pair<float, float> Pj = project2D(vertices[j]);

        // Skip degenerate edges (both vertices same Y)
        if (fabs(Pj.second - Pi.second) < 1e-8f)
            continue;

        // Check if the point is between y-bounds of the edge
        bool intersectY = ((Pi.second > P.second) != (Pj.second > P.second));
        if (intersectY) {
            float slope = (Pj.first - Pi.first) / (Pj.second - Pi.second);
            float xIntersect = slope * (P.second - Pi.second) + Pi.first;
            if (P.first < xIntersect)
                inside = !inside;
        }
    }

    return inside;
}


void Plane::sortVerticesWinding(){
    //Compute centroid of all vertices
    vector<float> center(3, 0.0f);
    for (const auto& v : vertices) {
        center[0] += v[0];
        center[1] += v[1];
        center[2] += v[2];
    }
    center[0] /= vertices.size();
    center[1] /= vertices.size();
    center[2] /= vertices.size();

    //Find the dominant plane axis
    int dropAxis;
    float ax = fabs(normal[0]);
    float ay = fabs(normal[1]);
    float az = fabs(normal[2]);
    if (ax > ay && ax > az) dropAxis = 0; // drop X
    else if (ay > az)       dropAxis = 1; // drop Y
    else                    dropAxis = 2; // drop Z

    //Compute 2D angles around the centroid
    struct VertexAngle {
        vector<float> vertex;
        float angle;
    };
    vector<VertexAngle> temp;
    temp.reserve(vertices.size());

    for (const auto& v : vertices) {
        float x, y, cx, cy;
        if (dropAxis == 0) { x = v[1]; y = v[2]; cx = center[1]; cy = center[2]; } // drop X
        else if (dropAxis == 1) { x = v[0]; y = v[2]; cx = center[0]; cy = center[2]; } // drop Y
        else { x = v[0]; y = v[1]; cx = center[0]; cy = center[1]; } // drop Z

        float angle = atan2(y - cy, x - cx);
        temp.push_back({v, angle});
    }

    //Sort by angle (CCW winding)
    sort(temp.begin(), temp.end(), [](const VertexAngle& a, const VertexAngle& b) {
        return a.angle < b.angle;
    });

    //Write back sorted vertices
    for (size_t i = 0; i < vertices.size(); ++i)
        vertices[i] = temp[i].vertex;
}


void Plane::calculateNormal(){
    float a[3] = {vertices[1][0] - vertices[0][0], vertices[1][1] - vertices[0][1], vertices[1][2] - vertices[0][2]};
    float b[3] = {vertices[2][0] - vertices[0][0], vertices[2][1] - vertices[0][1], vertices[2][2] - vertices[0][2]};

    normal[0] = a[1]*b[2] - a[2]*b[1]; 
    normal[1] = a[2]*b[0] - a[0]*b[2]; 
    normal[2] = a[0]*b[1] - a[1]*b[0];

    float normal_length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

    normal[0] /= normal_length;
    normal[1] /= normal_length;
    normal[2] /= normal_length;
}


AABB Plane::getAABB() const {
    vector<float> min = {FLT_MAX, FLT_MAX, FLT_MAX};
    vector<float> max = {FLT_MIN, FLT_MIN, FLT_MIN};

    // Calculate the minimum and maximum points of the plane's AABB
    for (int i = 0; i < vertices.size(); i++) {
        vector<float> vertex = vertices[i]; // Assuming vertices is a vector<vector<float>> member variable
        if (vertex[0] < min[0]) min[0] = vertex[0];
        if (vertex[0] > max[0]) max[0] = vertex[0];
        if (vertex[1] < min[1]) min[1] = vertex[1];
        if (vertex[1] > max[1]) max[1] = vertex[1];
        if (vertex[2] < min[2]) min[2] = vertex[2];
        if (vertex[2] > max[2]) max[2] = vertex[2];
    }

    return AABB(min, max);
}

