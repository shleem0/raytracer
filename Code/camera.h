#ifndef CAMERA
#define CAMERA

#include <string>
#include <vector>
#include "ray.h"

using namespace std;

class Camera{
    public:
        float location[3];
        float gaze_vector[3];
        float focal_length;
        float sensor_width;
        float sensor_height;
        float res_x;
        float res_y;

        static vector<Camera> parseCameraDataFromJson();
        Ray convertPixelToRay(float, float) const;
        static float getFloat(const string&, const string&);
        static string getJSONObject(const string&, const string&);
        static string getJSONArray(const string&, const string&);
};
#endif