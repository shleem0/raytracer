#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "camera.h"
#include "ray.h"
#include "shape.h"
#include "raytracer.h"

using namespace std;


//Overarching function for parsing camera data from JSON
vector<Camera> Camera::parseCameraDataFromJson(Config config) {
    //Open JSON file
    ifstream jsonFile("..\\ASCII\\scene.json");
    if (!jsonFile) {
        cerr << "Error opening JSON file." << endl;
        return {};
    }

    //Read entire file into string
    stringstream buffer;
    buffer << jsonFile.rdbuf();
    string jsonStr = buffer.str();

    //Get cameras array as string
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string camerasArrayStr = getJSONArray(propertiesStr, "\"cameras\"");

    //Remove outer brackets
    if (camerasArrayStr.front() == '[' && camerasArrayStr.back() == ']') {
        camerasArrayStr = camerasArrayStr.substr(1, camerasArrayStr.size() - 2);
    }

    //Split array into full camera objects
    vector<string> cameraObjects;
    int braces = 0;
    size_t start = 0;
    for (size_t i = 0; i < camerasArrayStr.size(); ++i) {
        if (camerasArrayStr[i] == '{') {
            if (braces == 0) start = i;
            braces++;
        } else if (camerasArrayStr[i] == '}') {
            braces--;
            if (braces == 0) {
                cameraObjects.push_back(camerasArrayStr.substr(start, i - start + 1));
            }
        }
    }

    if (cameraObjects.empty()){
        cerr << "No cameras." << endl;
        return {};
    }

    vector<Camera> cameras;


    for (int i = 0; i < cameraObjects.size(); i++){
        Camera newCam;

        //Extract the specific camera
        string cameraDataStr = cameraObjects[i];

        //Parse location
        string locationStr = getJSONObject(cameraDataStr, "\"location\"");
        newCam.location = {getFloat(locationStr, "\"x\""),
        getFloat(locationStr, "\"y\""),
        getFloat(locationStr, "\"z\"")};

        //Parse gaze vector
        string gazeVectorStr = getJSONObject(cameraDataStr, "\"gaze_vector\"");
        newCam.gaze_vector = {getFloat(gazeVectorStr, "\"x\""),
        getFloat(gazeVectorStr, "\"y\""),
        getFloat(gazeVectorStr, "\"z\"")};

        //Parse camera aperture
        if (config.dof){
            newCam.aperture = getFloat(cameraDataStr, "\"aperture\"");
        }
        else{
            newCam.aperture = 0.0f;
        }

        //Parse focal length
        newCam.focal_length = getFloat(cameraDataStr, "\"focal_length\"") / 1000.0f;
        newCam.focal_distance = getFloat(cameraDataStr, "\"focal_distance\"");

        //Parse sensor width and height
        string sensorStr = getJSONObject(cameraDataStr, "\"sensor\"");
        newCam.sensor_width = getFloat(sensorStr, "\"width\"") / 1000.0f;
        newCam.sensor_height = getFloat(sensorStr, "\"height\"") / 1000.0f;

        //Parse resolution
        string filmResolutionStr = getJSONObject(cameraDataStr, "\"film_resolution\"");
        newCam.res_x = getFloat(filmResolutionStr, "\"width\"");
        newCam.res_y = getFloat(filmResolutionStr, "\"height\"");

        cameras.push_back(newCam);
    }
    return cameras;
}


//Compute ray for pixel (x, y)
Ray Camera::convertPixelToRay(float x, float y) const {

    Ray ray;
    ray.time = static_cast<float>(rand()) / RAND_MAX;

    //Normalize pixel coordinates [0,1]
    float u = (x + 0.5f) / res_x;
    float v = (y + 0.5f) / res_y;

    //Convert to camera space
    float cam_x = -(u - 0.5f) * (sensor_width / focal_length);
    float cam_y = (0.5f - v) * (sensor_height / focal_length);
    float cam_z = -1.0f; //Blender camera looks along -Z

    //Normalize forward
    vector<float> forward = Raytracer::normalise(gaze_vector);

    //Compute right vector
    vector<float> world_up = {0.0f, 0.0f, 1.0f};

    // Standard right-handed camera: right = forward × world_up
    vector<float> right = {
        forward[1]*world_up[2] - forward[2]*world_up[1],
        forward[2]*world_up[0] - forward[0]*world_up[2],
        forward[0]*world_up[1] - forward[1]*world_up[0]
    };

    float rlen = sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);

    if (rlen < 1e-6f) {
        //forward was almost parallel to world_up, choose alternative up
        world_up[0] = 0.0f; world_up[1] = 0.0f; world_up[2] = 1.0f;
        right[0] = world_up[1]*forward[2] - world_up[2]*forward[1];
        right[1] = world_up[2]*forward[0] - world_up[0]*forward[2];
        right[2] = world_up[0]*forward[1] - world_up[1]*forward[0];
    }

    right = Raytracer::normalise(right);

    //Compute up vector
    vector<float> up = {
        forward[1]*right[2] - forward[2]*right[1],
        forward[2]*right[0] - forward[0]*right[2],
        forward[0]*right[1] - forward[1]*right[0]
    };

    //Compute direction in world space
    ray.direction = {cam_x*right[0] + cam_y*up[0] + cam_z*forward[0],
        cam_x*right[1] + cam_y*up[1] + cam_z*forward[1],
        cam_x*right[2] + cam_y*up[2] + cam_z*forward[2]};

    //Normalize
    float dlen = sqrt(ray.direction[0]*ray.direction[0] +
                        ray.direction[1]*ray.direction[1] +
                        ray.direction[2]*ray.direction[2]);
    if (dlen < 1e-6f) {
        //Handle division by zero or very small value, e.g., set a default direction
        ray.direction[0] = 0.0f; ray.direction[1] = 0.0f; ray.direction[2] = 1.0f;
        dlen = sqrt(ray.direction[0]*ray.direction[0] +
                        ray.direction[1]*ray.direction[1] +
                        ray.direction[2]*ray.direction[2]);
    }

    ray.direction = Raytracer::normalise(ray.direction);

    ray.direction[0] = -ray.direction[0];
    ray.direction[1] = -ray.direction[1];
    ray.direction[2] = -ray.direction[2];
    
    ray.origin = {location[0], location[1], location[2]};

    //Adjust ray for depth of field
    if (aperture > 0.0f){
        float lens_radius = focal_length/(2.0f * aperture);

        float sx, sy;
        sampleDisk(sx, sy);

        sx *= lens_radius;
        sy *= lens_radius;

        vector<float> lensOffset = Raytracer::add_vec(Raytracer::mul_vec(right, sx), Raytracer::mul_vec(up, sy));
        vector<float> focusPoint = Raytracer::add_vec(ray.origin, Raytracer::mul_vec(ray.direction, focal_distance));

        ray.origin = Raytracer::add_vec(ray.origin, lensOffset);

        ray.direction = Raytracer::normalise(Raytracer::sub_vec(focusPoint, ray.origin));
    }

    return ray;
}



/*Helper functions for parsing JSON data*/
//Function to get a float value from a JSON string
float Camera::getFloat(const string& str, const string& key) {
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find(':', pos);
    if (start == string::npos) {
        throw runtime_error("Invalid JSON string");
    }

    start += 1; //Skip the colon
    while (start < str.length() && isspace(str[start])) {
        start++;
    }

    size_t end = start;
    while (end < str.length() && isdigit(str[end]) || str[end] == '.' || str[end] == '-' || str[end] == 'e') {
        end++;
    }

    return stof(str.substr(start, end - start));
}

//Function to get a JSON object value from a JSON string
string Camera::getJSONObject(const string& str, const string& key) {
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find('{', pos);
    if (start == string::npos) {
        throw runtime_error("Invalid JSON string");
    }

    size_t end = start;
    int braces = 1;
    while (end < str.length()) {
        if (str[end] == '{') {
            braces++;
        } else if (str[end] == '}') {
            braces--;
            if (braces == 0) {
                break;
            }
        }
        end++;
    }

    return str.substr(start, end - start + 1);
}


//Function to get a JSON array from a JSON string
string Camera::getJSONArray(const string& str, const string& key) {
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find('[', pos);
    if (start == string::npos) {
        throw runtime_error("Invalid JSON string: '[' not found");
    }

    int brackets = 1;
    size_t end = start + 1;
    while (end < str.size() && brackets > 0) {
        if (str[end] == '[') brackets++;
        else if (str[end] == ']') brackets--;
        end++;
    }

    if (brackets != 0) {
        throw runtime_error("Mismatched brackets in JSON string");
    }

    return str.substr(start, end - start); //includes outer [ ... ]
}

void Camera::sampleDisk(float &sx, float &sy){


    // Map to [-1,1]
    float x = 2.0f * (static_cast<float>(rand()) / RAND_MAX) - 1.0f;
    float y = 2.0f * (static_cast<float>(rand()) / RAND_MAX) - 1.0f;

    if (x == 0.0f && y == 0.0f) {
        sx = 0.0f;
        sy = 0.0f;
        return;
    }

    float theta;
    float r;
    if (fabs(x) > fabs(y)){
        r = x;
        theta = (M_PI / 4.0f) * (y/x);
    }
    else{
        r = y;
        theta = (M_PI / 2.0f) - (M_PI / 4.0f) * (x / y);
    }

    sx = r * std::cos(theta);
    sy = r * std::sin(theta);
}