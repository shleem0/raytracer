#ifndef CONFIG_H
#define CONFIG_H

using namespace std;

struct Config {
    bool softShadows = false;
    int SSsamples = 4;
    bool glossyReflect = false;
    int GRsamples = 16;
    float lightRadius = 0.25f;

    bool antiAliasing = false;
    int AAsamples = 4;

    bool textures = false;

    bool bvh = true;
    
    bool reflections = false;
    int reflectRefractDepth = 0;

    bool dof = false;
    int dofSamples = 16;
    bool motionBlur = false;
    int mbSamples = 16;

    static Config parseArgs(int argc, char* argv[]);
};

#endif