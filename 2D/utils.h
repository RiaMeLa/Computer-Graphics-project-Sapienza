#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>


inline float clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

float lerp(const float &a, const float &b, const float &t)
{
    return a * (1 - t) + b * t;
}




inline float random_float() {
    return rand() / (RAND_MAX + 1.0);
}

inline float random_float(float min, float max) {
    return min + (max-min)*random_float();
}  

float fract(float x) {
    float n = std::floor(x);
    return x - n;
}

#endif