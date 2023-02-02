#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include "vec3.h"


float lerp(const float &a, const float &b, const float &t)
{
    return a * (1 - t) + b * t;
}

vec3 lerp(const vec3 &v,const vec3 &u, const float &t) {
    float x = lerp(v.x(), u.x(), t);
    float y = lerp(v.y(), u.y(), t);
    float z = lerp(v.z(), u.z(), t);
    return vec3(x, y, z);
}

/*inline float random_float() {
    return rand() / (RAND_MAX + 1.0);
}

inline float random_float(float min, float max) {
    return min + (max-min)*random_float();
}  */



#endif