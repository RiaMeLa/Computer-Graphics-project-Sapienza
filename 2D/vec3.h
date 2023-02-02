#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "utils.h"

inline float random_float(float min, float max);

class vec3 {
    public:
        vec3() : x(0), y(0), z(0) {}
        vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
        vec3& operator *= (const float &r) { x *= r, y *= r, z *= r; return *this; }
        vec3& operator += (const float &r) { x += r, y += r, z += r; return *this; }
        vec3& operator += (const vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
        vec3& operator /= (const float &r) { return *this *= float(1/r); }
        float length_squared() { return x * x + y * y + z * z; }
        float length() { return sqrt(length_squared());}
        vec3& normalize() { *this /= length(); return *this; }

        inline static vec3 random(float min, float max) { 
            return vec3(random_float(min,max), random_float(min,max), random_float(min,max));
        }
    public:
        float x;
        float y;
        float z;
};

inline std::ostream& operator << (std::ostream &out, vec3 &v) {
    return out << v.x << ' ' << v.y << ' ' << v.z;
}

inline vec3 operator + (vec3 &u, vec3 &v) { 
    return vec3(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline vec3 operator - (vec3 &u, vec3 &v) { 
    return vec3(u.x - v.x, u.y - v.y, u.z - v.z);
}

inline vec3 operator - (vec3 u, float s) { 
    return vec3(u.x - s, u.y - s, u.z - s);
}

inline vec3 operator * (vec3 u, vec3 v) { 
    return vec3(u.x * v.x, u.y * v.y, u.z * v.z);
}

inline vec3 operator * (float s, vec3 v) {
    return vec3(s*v.x, s*v.y, s*v.z);
}

inline vec3 operator * (const vec3 &v, const float &s) {
    return s*v;
}

inline vec3 operator + (float s, vec3 v) {
    return vec3(s*v.x, s*v.y, s*v.z);
}

inline vec3 operator/ (const vec3 &v, float s) {
    return (1/s) * v;
}

inline float dot(vec3 &u, vec3 &v) {
    return u.x * v.x
         + u.y * v.y
         + u.z * v.z;
}

inline vec3 cross(vec3 &u, vec3 &v) {
     return vec3(u.y * v.z - u.z * v.y,
                 u.z * v.x - u.x * v.z,
                 u.x * v.y - u.y * v.x);
}

vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1,1);
        if (p.length() >= 1) continue;
        return p;
    }
}

vec3 random_unit_vec3() {
    vec3 v = random_in_unit_sphere();
    v.normalize();
    return v;
}

vec3 sin(vec3 v) {
    vec3 u(sinf(v.x), sinf(v.y), sinf(v.z));
    return u;
}

vec3 sqrt(vec3 v) {
    return vec3(sqrt(v.x), sqrt(v.y), sqrt(v.z));
}

vec3 lerp(const vec3 &v,const vec3 &u, const float &t) {
    float x = lerp(v.x, u.x, t);
    float y = lerp(v.y, u.y, t);
    float z = lerp(v.z, u.z, t);
    return vec3(x, y, z);
}
#endif