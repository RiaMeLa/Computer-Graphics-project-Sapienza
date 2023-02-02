#ifndef VEC2_H
#define VEC2_H

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "utils.h"

class vec2 {
    public:
        vec2() : x(0), y(0) {}
        vec2(float xx, float yy) : x(xx), y(yy) {}
        vec2& operator *= (float r) { x *= r, y *= r; return *this; }

    public:
        float x;
        float y;
};


inline vec2 operator + (vec2 u, vec2 v) { 
    return vec2(u.x + v.x, u.y + v.y);
}

inline vec2 operator + (float n, vec2 v) { 
    return vec2(v.x + n, v.y + n);
}

inline vec2 operator + (vec2 v, float n) { 
    return vec2(v.x + n, v.y + n);
}

inline vec2 operator - (vec2 u, vec2 v) { 
    return vec2(u.x - v.x, u.y - v.y);
}

inline vec2 operator - (vec2 u, float s) { 
    return vec2(u.x - s, u.y - s);
}

inline vec2 operator * (const vec2 &u,const vec2 &v) { 
    return vec2(u.x * v.x, u.y * v.y);
}

inline vec2 operator * (float s, vec2 v) {
    return vec2(s*v.x, s*v.y);
}

inline vec2 operator * (vec2 v, float s) {
    return s*v;
}

inline vec2 operator / (vec2 v, float s) {
    return (1/s) * v;
}

inline float dot(vec2 u, vec2 v) {
    return u.x * v.x
         + u.y * v.y;
}

vec2 fract(vec2 v) {
    vec2 u(fract(v.x), fract(v.y));
    return u;
}

vec2 sin(vec2 v) {
    return vec2(sinf(v.x), sinf(v.y));
}

vec2 floor(vec2 v) {
    return vec2(std::floor(v.x), std::floor(v.y));
}

vec2 normalize(vec2 v) {
    float length = sqrt(v.x * v.x + v.y * v.y);
    return v / length;
}

vec2 lerp(const vec2 &v,const vec2 &u, const float &t) {
    float x = lerp(v.x, u.x, t);
    float y = lerp(v.y, u.y, t);
    return vec2(x, y);
}


#endif