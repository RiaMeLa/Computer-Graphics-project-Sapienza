#ifndef REMAPPERS_H
#define REMAPPERS_H


inline float smoothstep(float t) {
    return t * t * (3 - 2 * t);
}

inline float smoothstep(float edge0, float edge1, float t) {
    t = clamp((t - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
    return t * t * (3 - 2 * t);
}

inline float quintic(float t) {
    return (( 6.0f * t - 15.0f ) * t + 10.0f) * t * t * t;
}


#endif