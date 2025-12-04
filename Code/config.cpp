#include <string>
#include <iostream>
#include "config.h"

using namespace std;

Config Config::parseArgs(int argc, char* argv[]){

    Config config;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        //Soft shadows
        if (arg == "-ss" || arg == "--soft_shadows") {
            config.softShadows = true;
        }
        else if ((arg == "-sss" || arg == "--ss_samples") && i + 1 < argc){
            config.SSsamples = std::stoi(argv[++i]);
        } 


        //Glossy reflections
        else if (arg == "-gr" || arg == "--glossy_reflect"){
            config.glossyReflect = true;
        }
        else if((arg == "-grs" || arg == "--gr_samples") && i + 1 < argc){
            config.GRsamples = std::stoi(argv[++i]);
        }


        //Anti-aliasing
        else if (arg == "-aa" || arg == "--antialiasing") {
            config.antiAliasing = true;
        }
        else if ((arg == "-aas" || arg == "--aa_samples") && i + 1 < argc) {
            config.AAsamples = std::stoi(argv[++i]);
        }


        //Deactivating BVH acceleration
        else if (arg == "-u" || arg == "--unaccelerated"){
            config.bvh = false;
        }


        //Reflections/refractions + depth
        else if (arg == "-r" || arg == "--reflections"){
            config.reflectRefractDepth = 1;
        }
        else if ((arg == "-rd" || arg == "--reflect_depth") && i + 1 < argc){
            config.reflectRefractDepth = std::stoi(argv[++i]);
        }


        //Texture mapping
        else if (arg == "-t" || arg == "--texture_mapping"){
            config.textures = true;
        }


        //Depth of field + samples
        else if (arg == "-dof" || arg == "--depthoffield"){
            config.dof = true;
        }

        else if ((arg == "-dofs" || arg == "--dof_samples") && i + 1 < argc){
            config.dofSamples = std::stoi(argv[++i]);
        }

        //Motion blur + samples
        else if ((arg == "-m" || arg == "--motion-blur")){
            config.motionBlur = true;
        }

        else if ((arg == "-mbs" || arg == "--mb_samples") && i + 1 < argc){
            config.mbSamples = std::stoi(argv[++i]);
        }

        //Output file name
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc){
            config.outputFile = argv[++i];
        }

        else {
            std::cerr << "Unknown flag: " << arg << "\n";
        }
    }

    return config;
}