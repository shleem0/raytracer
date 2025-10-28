#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include "cube.h"
#include "ray.h"

using namespace std;

static void Transformer::rotateRayInverse(const Ray& ray, Ray& localRay, Cube *cube){
    // Inverse rotation (reverse order, negate angles)
    float cx = cos(-rot[0]), sx = sin(-rot[0]);
    float cy = cos(-rot[1]), sy = sin(-rot[1]);
    float cz = cos(-rot[2]), sz = sin(-rot[2]);

    // Copy origin and direction
    out = ray;

    auto rotate = [&](array<float,3>& v) {
        float x = v[0], y = v[1], z = v[2];

        // Rotate around Z
        float xz = x*cz - y*sz;
        float yz = x*sz + y*cz;
        x = xz; y = yz;

        // Rotate around Y
        float xy = x*cy + z*sy;
        float zy = -x*sy + z*cy;
        x = xy; z = zy;

        // Rotate around X
        float yy = y*cx - z*sx;
        float zy2 = y*sx + z*cx;
        y = yy; z = zy2;

        v[0] = x; v[1] = y; v[2] = z;
    };

    // Translate origin into cube space (inverse translation)
    out.origin[0] -= ray.origin[0];
    out.origin[1] -= ray.origin[1];
    out.origin[2] -= ray.origin[2];

    rotate(out.origin);
    rotate(out.direction);
}



static void transRayToLocal(const Ray& ray, Ray& local, Cube *cube){
    // Move origin relative to cube
    local.origin[0] = (ray.origin[0] - location[0]) / cube.scale;
    local.origin[1] = (ray.origin[1] - location[1]) / cube.scale;
    local.origin[2] = (ray.origin[2] - location[2]) / cube.scale;

    local.direction = ray.direction;
    // Apply inverse rotation
    float rotInv[3] = {-cube.rotation[0], -cube.rotation[1], -cube.rotation[2]};
    Ray tmp = local;
    rotateRayInverse(tmp, local, rotInv);

    // Scale direction too
    local.direction[0] /= scale;
    local.direction[1] /= scale;
    local.direction[2] /= scale;
}

static bool intersectUnitCube(const Ray& ray, float tNear, float tFar) {
    float tmin = -INFINITY, tmax = INFINITY;

    for (int i = 0; i < 3; ++i) {
        if (fabs(ray.direction[i]) < 1e-6) {
            // Ray is parallel to slab
            if (ray.origin[i] < -1 || ray.origin[i] > 1)
                return false;
        } else {
            float t1 = (-1 - ray.origin[i]) / ray.direction[i];
            float t2 = ( 1 - ray.origin[i]) / ray.direction[i];
            if (t1 > t2) swap(t1, t2);
            tmin = max(tmin, t1);
            tmax = min(tmax, t2);
            if (tmin > tmax)
                return false;
        }
    }
    tNear = tmin;
    tFar = tmax;
    return true;
};