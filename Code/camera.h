#ifndef CAMERA
#define CAMERA

#include <string>
#include <vector>
#include "ray.h"

using namespace std;

//Class for a Blender camera object
class Camera{
    public:
        vector<float> location;
        vector<float> gaze_vector;
        float focal_length;
        float sensor_width;
        float sensor_height;
        float res_x;
        float res_y;

        //Parses camera data from the scene's JSON into the object
        static vector<Camera> parseCameraDataFromJson();

        //Converts a given pixel (x,y) to a ray
        Ray convertPixelToRay(float, float) const;

        //Extract a float value from a JSON object
        static float getFloat(const string&, const string&);

        //Extract a value from a JSON object
        static string getJSONObject(const string&, const string&);

        //Extract an array from a JSON object
        static string getJSONArray(const string&, const string&);
};
#endif