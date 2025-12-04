To run the raytracer:
- Export a Blender scene by running Export.py within Blender
- Compile code using the makefile by running "make"
- Run the raytracer using "raytracer.exe"

To customise raytracing, run with these commands after raytracer.exe:

- "-s" or "--soft_shadows" to enable soft shadows
- "-sss" or "-ss_samples" followed by a number to pick number of soft shadow samples (default 4)

- "-gr" or "--glossy_reflect" to enable glossy reflections (must be enabled with regular reflections)
- "-grs" or "--gr_samples" followed by a number to pick number of glossy reflection samples (default 16)

- "-aa" or "--antialiasing" to enable anti-aliasing
- "-aas" or "--aa_samples" followed by a number to pick number of anti-aliasing samples (default 4)

- "-u" or "--unaccelerated" to disable BVH-accelerated raytracing

- "-r" or "--reflections" to enable reflection and refraction
- "-rd" or "--reflect_depth" followed by a number to pick reflection/refraction depth (default 1)

- "-t" or "--texture_mapping" to enable texture mapping

- "-dofs" or "--depthoffield" to enable depth of field
- "-dofs" or "--dof_samples" followed by a number to pick number of depth of field samples (default 16)

- "-m" or "--motion_blur" to enable motion blur
- "-mbs" or "--mb_samples" followed by a number to pick the number of motion blur samples (default 16)

- "-o" or "--output" followed by a string (ending in .ppm) to decide output image file name