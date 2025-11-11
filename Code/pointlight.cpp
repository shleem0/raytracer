#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include "pointlight.h"


using namespace std;

vector<PointLight> PointLight::parseLightDataFromJson() {
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

    //Get lights array as string
    string propertiesStr = getJSONObject(jsonStr, "\"properties\"");
    string lightArrayStr = getJSONArray(propertiesStr, "\"point_lights\"");

    //Remove outer brackets
    if (lightArrayStr.front() == '[' && lightArrayStr.back() == ']') {
        lightArrayStr = lightArrayStr.substr(1, lightArrayStr.size() - 2);
    }

    //Split array into full light objects
    vector<string> lightObjects;
    int braces = 0;
    size_t start = 0;
    for (size_t i = 0; i < lightArrayStr.size(); ++i) {
        if (lightArrayStr[i] == '{') {
            if (braces == 0) start = i;
            braces++;
        } else if (lightArrayStr[i] == '}') {
            braces--;
            if (braces == 0) {
                lightObjects.push_back(lightArrayStr.substr(start, i - start + 1));
            }
        }
    }

    if (lightObjects.empty()){
        cerr << "No lights." << endl;
        return {};
    }

    vector<PointLight> lights;



    for (int i = 0; i < lightObjects.size(); i++){
        PointLight newLight;

        //Extract the specific light
        string lightDataStr = lightObjects[i];

        //Parse location
        string locationStr = getJSONObject(lightDataStr, "\"location\"");
        newLight.location[0] = getFloat(locationStr, "\"x\"");
        newLight.location[1] = getFloat(locationStr, "\"y\"");
        newLight.location[2] = getFloat(locationStr, "\"z\"");

        //Parse radiant intensity
        newLight.rad_intensity = getFloat(lightDataStr, "\"radiant_intensity\"");

        lights.push_back(newLight);
    }
    return lights;
}

/*Helper functions for parsing JSON data*/
//Function to get a float value from a JSON string
float PointLight::getFloat(const string& str, const string& key) {
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
    while (end < str.length() && isdigit(str[end]) || str[end] == '.' || str[end] == '-') {
        end++;
    }

    return stof(str.substr(start, end - start));
}

//Function to get a JSON object value from a JSON string
string PointLight::getJSONObject(const string& str, const string& key) {
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
string PointLight::getJSONArray(const string& str, const string& key) {
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