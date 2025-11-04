#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include "shape.h"

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
    while (end < str.length() && isdigit(str[end]) || str[end] == '.' || str[end] == '-') {
        end++;
    }

    return stof(str.substr(start, end - start));
}

//Function to get a JSON object value from a JSON string
string Shape::getJSONObject(const string& str, const string& key) {
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
