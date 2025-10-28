#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "camera.h"
#include "ray.h"
#include "shape.h"

using namespace std;


//Overarching function for parsing camera data from JSON
vector<Camera> Camera::parseCameraDataFromJson() {
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

    // Get cameras array as string (includes [ ... ])
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string camerasArrayStr = getJSONArray(propertiesStr, "\"cameras\"");

    // Remove outer brackets
    if (camerasArrayStr.front() == '[' && camerasArrayStr.back() == ']') {
        camerasArrayStr = camerasArrayStr.substr(1, camerasArrayStr.size() - 2);
    }

    // Split array into full camera objects
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

        // Extract the specific camera
        string cameraDataStr = cameraObjects[i];

        // Parse location
        string locationStr = getJSONObject(cameraDataStr, "\"location\"");
        newCam.location[0] = getFloat(locationStr, "\"x\"");
        newCam.location[1] = getFloat(locationStr, "\"y\"");
        newCam.location[2] = getFloat(locationStr, "\"z\"");

        // Parse gaze vector
        string gazeVectorStr = getJSONObject(cameraDataStr, "\"gaze_vector\"");
        newCam.gaze_vector[0] = getFloat(gazeVectorStr, "\"x\"");
        newCam.gaze_vector[1] = getFloat(gazeVectorStr, "\"y\"");
        newCam.gaze_vector[2] = getFloat(gazeVectorStr, "\"z\"");

        // Parse focal length
        newCam.focal_length = getFloat(cameraDataStr, "\"focal_length\"") / 1000.0f;

        // Parse sensor size
        string sensorStr = getJSONObject(cameraDataStr, "\"sensor\"");
        newCam.sensor_width = getFloat(sensorStr, "\"width\"") / 1000.0f;
        newCam.sensor_height = getFloat(sensorStr, "\"height\"") / 1000.0f;

        // Parse film resolution
        string filmResolutionStr = getJSONObject(cameraDataStr, "\"film_resolution\"");
        newCam.res_x = getFloat(filmResolutionStr, "\"width\"");
        newCam.res_y = getFloat(filmResolutionStr, "\"height\"");

        cameras.push_back(newCam);
    }
    return cameras;
}


// Compute ray for pixel (x, y)
Ray Camera::convertPixelToRay(float x, float y) const {

    Ray ray;

    //Normalize pixel coordinates [0,1]
    float u = (x + 0.5f) / res_x;
    float v = (y + 0.5f) / res_y;

    //Convert to camera space
    float cam_x = (u - 0.5f) * (sensor_width / focal_length);
    float cam_y = (0.5f - v) * (sensor_height / focal_length);
    float cam_z = -1.0f; // Blender camera looks along -Z

    //Normalize forward
    float len = sqrt(gaze_vector[0]*gaze_vector[0] +
                        gaze_vector[1]*gaze_vector[1] +
                        gaze_vector[2]*gaze_vector[2]);

    float forward[3] = { gaze_vector[0]/len, gaze_vector[1]/len, gaze_vector[2]/len };

    //Compute right vector
    //cross(world_up, forward) to ensure right-hand system
    float world_up[3] = {0.0f, 1.0f, 0.0f};
    float right[3] = {
        world_up[1]*forward[2] - world_up[2]*forward[1],
        world_up[2]*forward[0] - world_up[0]*forward[2],
        world_up[0]*forward[1] - world_up[1]*forward[0]
    };
    float rlen = sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);

    if (rlen < 1e-6f) {
        // forward was almost parallel to world_up, choose alternative up
        world_up[0] = 0.0f; world_up[1] = 0.0f; world_up[2] = 1.0f;
        right[0] = world_up[1]*forward[2] - world_up[2]*forward[1];
        right[1] = world_up[2]*forward[0] - world_up[0]*forward[2];
        right[2] = world_up[0]*forward[1] - world_up[1]*forward[0];
        rlen = sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);
    }
    right[0] /= rlen; right[1] /= rlen; right[2] /= rlen;

    //Compute up vector
    float up[3] = {
        forward[1]*right[2] - forward[2]*right[1],
        forward[2]*right[0] - forward[0]*right[2],
        forward[0]*right[1] - forward[1]*right[0]
    };

    //Compute direction in world space
    ray.direction[0] = cam_x*right[0] + cam_y*up[0] + cam_z*forward[0];
    ray.direction[1] = cam_x*right[1] + cam_y*up[1] + cam_z*forward[1];
    ray.direction[2] = cam_x*right[2] + cam_y*up[2] + cam_z*forward[2];

    //Normalize
    float dlen = sqrt(ray.direction[0]*ray.direction[0] +
                        ray.direction[1]*ray.direction[1] +
                        ray.direction[2]*ray.direction[2]);
    if (dlen < 1e-6f) {
        // Handle division by zero or very small value, e.g., set a default direction
        ray.direction[0] = 0.0f; ray.direction[1] = 0.0f; ray.direction[2] = 1.0f;
        dlen = sqrt(ray.direction[0]*ray.direction[0] +
                        ray.direction[1]*ray.direction[1] +
                        ray.direction[2]*ray.direction[2]);
    }
    ray.direction[0] /= dlen;
    ray.direction[1] /= dlen;
    ray.direction[2] /= dlen;

    ray.direction[0] = -ray.direction[0];
    ray.direction[1] = -ray.direction[1];
    ray.direction[2] = -ray.direction[2];
    

    ray.origin[0] = location[0];
    ray.origin[1] = location[1];
    ray.origin[2] = location[2];

    return ray;
}



/*Helper functions for parsing JSON data*/
// Function to get a float value from a JSON string
float Camera::getFloat(const string& str, const string& key) {
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find(':', pos);
    if (start == string::npos) {
        throw runtime_error("Invalid JSON string");
    }

    start += 1; // Skip the colon
    while (start < str.length() && isspace(str[start])) {
        start++;
    }

    size_t end = start;
    while (end < str.length() && isdigit(str[end]) || str[end] == '.' || str[end] == '-') {
        end++;
    }

    return stof(str.substr(start, end - start));
}

// Function to get a JSON object value from a JSON string
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

    return str.substr(start, end - start); // includes outer [ ... ]
}