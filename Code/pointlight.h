#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H

#include <vector>
#include <string>

using namespace std;

class PointLight{
    public:
        vector<float> location;
        float rad_intensity;

        static vector<PointLight> parseLightDataFromJson();
        static float getFloat(const string&, const string&);
        static string getJSONObject(const string&, const string&);
        static string getJSONArray(const string&, const string&);

};

#endif