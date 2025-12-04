#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "shape.h"
#include "raytracer.h"

using namespace std;

/*Helper functions for parsing JSON data*/
//Function to get a float value from a JSON string
float Shape::getFloat(const string& str, const string& key) {
    
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


 string Shape::getString(const string& str, const string& key) {
    
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find(':', pos);
    if (start == string::npos) {
        throw runtime_error("Invalid JSON string");
    }

    // Move just past colon and skip whitespace
    start = str.find_first_not_of(" \t\n\r", start + 1);
    if (start == string::npos) {
        return "";
    }

    // Case 1: null value
    if (str.compare(start, 4, "null") == 0) {
        return "";
    }

    // Case 2: quoted string value
    if (str[start] == '"') {
        size_t end = str.find('"', start + 1);
        if (end == string::npos) {
            throw runtime_error("Unterminated string value in JSON");
        }
        return str.substr(start + 1, end - start - 1);
    }

    // Otherwise unsupported format → return empty
    return "";
}

//Function to get a JSON object value from a JSON string
string Shape::getJSONObject(const string& str, const string& key) {
    size_t pos = str.find(key);
    if (pos == string::npos) {
        throw runtime_error("Key not found in JSON string");
    }

    size_t start = str.find('{', pos);
    if (start == string::npos) {
        if (key == "\"texture\""){
            return "";
        }
        else{
            throw runtime_error("Invalid JSON string");
        }
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

//Function to get an array from a JSON string
string Shape::getJSONArray(const string& str, const string& key) {
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

    return str.substr(start, end - start);
}

//Get linear interpolation of shape position for motion blur
vector<float> Shape::positionAt(float time){
    return {
        Raytracer::lerp(startLocation[0], endLocation[0], time),
        Raytracer::lerp(startLocation[1], endLocation[1], time),
        Raytracer::lerp(startLocation[2], endLocation[2], time)
    };
}
