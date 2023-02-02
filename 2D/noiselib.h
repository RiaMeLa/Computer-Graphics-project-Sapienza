#ifndef NOISELIB_H
#define NOISELIB_H

#include <stdlib.h>
#include <cstdio> 
#include <random> 
#include <functional> 
#include <iostream> 
#include <fstream> 
#include <cmath> 

#include "vec2.h"
#include "vec3.h"
#include "utils.h"
#include "remappers.h"



class value_noise {
    public:
        value_noise(int s) {
            seed = s;
            srand(seed);

            for (int k = 0; k < max_table_size; ++k){
                r[k] = random_float();
                permutation_table[k] = k;
            }

            for (int k = 0; k < max_table_size; ++k) {
                int i = rand() & max_size_mask;
                std::swap(permutation_table[k], permutation_table[i]);
                permutation_table[k + max_table_size] = permutation_table[k];
            }
        }

        float eval(vec2 p) {
            int xi = std::floor(p.x);
            int yi = std::floor(p.y);

            float tx = p.x - xi;
            float ty = p.y - yi;

            int xi0 = xi & max_size_mask;
            int xi1 = (xi0 + 1) & max_size_mask;
            int yi0 = yi & max_size_mask;
            int yi1 = (yi0 + 1) & max_size_mask;

            float c00 = r[permutation_table[permutation_table[xi0] + yi0]];
            float c01 = r[permutation_table[permutation_table[xi0] + yi1]];
            float c10 = r[permutation_table[permutation_table[xi1] + yi0]];
            float c11 = r[permutation_table[permutation_table[xi1] + yi1]];

        // remapping of tx and ty using the Smoothstep function 
            float sx = smoothstep(tx); 
            float sy = smoothstep(ty); 
    
            // linearly interpolate values along the x axis
            float nx0 = lerp(c00, c10, sx); 
            float nx1 = lerp(c01, c11, sx);

            return lerp(nx0, nx1, sy);
        }

    private:
        int seed;
        static const int max_table_size = 256;
        static const int max_size_mask  = max_table_size -1;
        float r[max_table_size];
        int permutation_table[max_table_size * 2];

    
};

class perlin_noise {
    public:
        perlin_noise (int s) {
            seed = s;
            srand(seed);
            for (int i = 0; i < table_size; ++i) { 
                gradients[i] = vec2(cosf(2.0f * M_PI * i / table_size), sin(2.0f * M_PI * i / table_size)); 
                permutation_table[i] = i; 
            }
            //shuffle permutation table
            for (int i = 0; i < table_size; ++i) 
                std::swap(permutation_table[i], permutation_table[rand() & mask_size]); 
            // extend the permutation table in the index range [256:512]
            for (int i = 0; i < table_size; ++i) 
                permutation_table[table_size + i] = permutation_table[i]; 
        }

        float eval(vec2 p) {  
            /*floor returns lowest int if input is signed*/
            int xi  = std::floor(p.x);
            int yi  = std::floor(p.y);

            int xi0 = xi & mask_size;
            int yi0 = yi & mask_size;
            int xi1 = (xi0 + 1) & mask_size;
            int yi1 = (yi0 + 1) & mask_size; 

            float tx = p.x - xi;
            float ty = p.y - yi;

            float sx = quintic(tx);
            float sy = quintic(ty);

            vec2 c00 = gradients[permutation_table[permutation_table[xi0] + yi0]];
            vec2 c01 = gradients[permutation_table[permutation_table[xi0] + yi1]];            
            vec2 c10 = gradients[permutation_table[permutation_table[xi1] + yi0]];
            vec2 c11 = gradients[permutation_table[permutation_table[xi1] + yi1]]; 
            
            //vectors from grid points to p
            float x0 = tx;
            float y0 = ty;
            float x1 = tx - 1;
            float y1 = ty - 1;

            vec2 p00 = vec2(x0,y0);           
            vec2 p01 = vec2(x0,y1);
            vec2 p10 = vec2(x1,y0);
            vec2 p11 = vec2(x1,y1);

            //lerp
            float a = lerp(dot(c00, p00), dot(c10, p10), sx);
            float b = lerp(dot(c01, p01), dot(c11, p11), sx);

            return (lerp(a, b, sy) + 1) * 0.5;
        }
    private:
        int seed;
        static const int table_size = 256; 
        static const int mask_size = table_size - 1; 
        vec2 gradients[table_size]; 
        unsigned permutation_table[table_size * 2]; 
};

class voronoi {
    public:
        voronoi() {}

        float eval(vec2 p) {
            vec2 n = vec2(std::floor(p.x), std::floor(p.y));
            vec2 f = fract(p);

            float res = 80.0;
            for (int j = -1; j <= 1; ++j)
                for (int i = -1; i <= 1; ++i) {
                    vec2 b = vec2(float(i), float(j));
                    //vec2 k = hash(b + n);
                    vec2 k = getPoint(n + b);
                    vec2 r = b - f + k;
                    //r = r + k;
                    float d = dot(r, r);
                    res = std::min(res, d);
                }
            return sqrt(res);
        }

    private:
        // Hash from "Hash without Sine" by Dave_Hoskins (https://www.shadertoy.com/view/4djSRW)
        float Hash11(float x) {
            x = fract(x * 0.1031);
            x *= x + 33.33;
            x *= x + x;
            return fract(x);
        }

        vec2 getPoint(vec2 cell) {
            float freq = Hash11(dot(cell, vec2(393.84, 673.48))) * 3.0 + 1.0;
            float phase = Hash11(dot(cell, vec2(348.46, 183.37)));
            float amp = Hash11(dot(cell, vec2(275.35, 741.69)));

            float t = freq + phase;
            return 0.5 + 0.5 * vec2(cos(t), sin(t)) * amp;
        }
};

class voronoi_colored {
    public:
        voronoi_colored() {}

        vec3 eval(vec2 p) {
            vec2 n = vec2(std::floor(p.x), std::floor(p.y));
            vec2 f = fract(p);

            float res = 8.0;
            vec3 col;
            for (int j = -1; j <= 1; ++j)
                for (int i = -1; i <= 1; ++i) {
                    vec2 b = vec2(float(i), float(j));
                    vec2 k = hash(b + n);
                    vec2 r = b - f + k;
                    //r = r + k;
                    float d = dot(r, r);

                    /*salvare i val di k, d scelti e computare fuori dal for*/
                    if (d < res) {
                        res = d;                         //k should be b + n
                        vec3 v = 0.5 + 0.5 * ( hash1(dot(k,vec2(7.0,113.0)))*1.5 + 1.8 + vec3(k.x + 2.7, k.x +  1.9, k.x + 1.7));
                        col = vec3(sinf(v.x) - d * 0.5,d * 0.75 - sinf(v.y), d * 0.75 - sinf(v.z));
                        col = col * col;
                    }
                }
            /*col = 0.5f + 0.5f * vec3(cos(col.x * 6.2831 + 0.0),
                                         cos(col.x * 6.2831 + 1.0),
                                         cos(col.x * 6.2831 + 2.0));*/

            return col;
        }

    private:
        vec2 hash(vec2 p) {
        //p = mod(p, 4.0); // tile
        p = vec2(dot(p,vec2(127.1,311.7)),
                dot(p,vec2(269.5,183.3)));
        p.x = sinf(p.x);
        p.y = sinf(p.y);
        return fract(p * 18.5453);
        /*int x = ((int) std::floor(p.x));
        int y = ((int) std::floor(p.y));
        return offsets[permutation_table[(int)(((int) x * 13.456) - ((int)y * 13.654)) & (mask_size * 2)]]*/;//[permutation_table[permutation_table[((int)p.x & mask_size)] + ((int) p.y & mask_size)]];
        }
        float hash1( float n ) { return fract(sin(n)*43758.5453); }
};

class smooth_voronoi {
    public:
        smooth_voronoi() {}

        float eval(vec2 p) {
            vec2 n = vec2(std::floor(p.x), std::floor(p.y));
            vec2 f = fract(p);

            float res = 0.0;
            for (int j = -1; j <= 1; ++j)
                for (int i = -1; i <= 1; ++i) {
                    vec2 b = vec2(float(i), float(j));
                    //vec2 k = hash(b + n);
                    vec2 k = getPoint(b + n);
                    vec2 r = b - f + k;
                    //r = r + k;
                    float d = dot(r, r);
                    res += 1.0/pow( d, 12.0 );
                }
            return pow( 1.0/res, 1.0/16.0 );
        }

    private:
        float Hash11(float x) {
            x = fract(x * 0.1031);
            x *= x + 33.33;
            x *= x + x;
            return fract(x);
        }

        vec2 getPoint(vec2 cell) {
            float freq = Hash11(dot(cell, vec2(393.84, 673.48))) * 3.0 + 1.0;
            float phase = Hash11(dot(cell, vec2(348.46, 183.37)));
            float amp = Hash11(dot(cell, vec2(275.35, 741.69)));

            float t = freq + phase;
            return 0.5 + 0.5 * vec2(cos(t), sin(t)) * amp;
        }

};

class voronoise {
    public: 
        voronoise() {}

        float eval(vec2 p, float u, float v) {
            vec2 n = vec2(std::floor(p.x), std::floor(p.y));
            vec2 f = fract(p);

            float sharpness = 1.0 + 63.0 * pow(1.0-v, 4.0);

            float value = 0.0;
            float accum = 0.0;
            for (int j = -2; j <= 2; ++j)
                for (int i = -2; i <= 2; ++i) {
                    vec2 curr = vec2(float(i), float(j));
                    vec2 cen = getPoint(curr + n, u);
                    vec2 r = curr - f + cen;
                    //vec2 k = hash(b + n);
                    //r = r + k;
                    float d = dot(r, r);
                    float w = pow( 1.0-smoothstep(0.0,1.414,sqrt(d)), sharpness);
                    float col = rand(n + curr);
                    value += w*col;
                    accum += w;
                }
            return value/accum;
        }

    private:
        float rand(vec2 p) {
            return fract(sin(dot(p,vec2(419.2,371.9))) * 833458.57832);
        }

        float Hash11(float x) {
            x = fract(x * 0.1031);
            x *= x + 33.33;
            x *= x + x;
            return fract(x);
        }

        vec2 getPoint(vec2 cell, float u) {
            float freq = Hash11(dot(cell, vec2(393.84, 673.48))) * 3.0 + 1.0;
            float phase = Hash11(dot(cell, vec2(348.46, 183.37)));
            float amp = Hash11(dot(cell, vec2(275.35, 741.69)));

            float t = freq + phase;
            return 0.5 + 0.5 * vec2(cos(t), sin(t)) * amp * u;
        }

};

class voronoi_distances {
    public:
        vec3 eval(vec2 p) {
            /*As first we search which cell contains the closest center 
            to our point p*/
            vec2 n = floor(p);
            vec2 f = fract(p);

            vec2 mg, mr;
            float md = 8.0;
            for(int j=-1; j<=1; j++)
            for(int i=-1; i<=1; i++) {
                vec2 curr = vec2(float(i), float(j));
                vec2 cen = getPoint(curr + n);
                vec2 r = curr + cen - f;
                float d = dot(r, r);

                if (d < md) {
                    md = d;
                    mr = r;
                    mg = curr;
                }
            }

            /*Next we search for the closest neighbor point around
             the cell we determined earlier*/
            md = 8.0;
            for(int j=-1; j<=1; j++)
            for(int i=-1; i<=1; i++) {
                vec2 g = mg + vec2(float(i), float(j));
                vec2 cen = getPoint(n + g);
                vec2 r = g + cen - f;
                if( dot(mr-r,mr-r)>0.00001 )
                    md = std::min(md, dot(0.5 * (mr + r), normalize(r - mr)));
            
            }
            return vec3(md, mr.x, mr.y);
        }

        float Hash11(float x) {
            x = fract(x * 0.1031);
            x *= x + 33.33;
            x *= x + x;
            return fract(x);
        }

        vec2 getPoint(vec2 cell) {
            float freq = Hash11(dot(cell, vec2(393.84, 673.48))) * 3.0 + 1.0;
            float phase = Hash11(dot(cell, vec2(348.46, 183.37)));
            float amp = Hash11(dot(cell, vec2(275.35, 741.69)));

            float t = freq + phase;
            return 0.5 + 0.5 * vec2(cos(t), sin(t)) * amp;
        }
};

class cloudy {
    public:

        float warp(vec2 p, float mm) {
            float m = 4.0;
            vec2 q = vec2(fbm(vec2(p)), fbm(p+vec2(5.12*0.01, 1.08)));
            
            vec2 r = vec2(fbm((p+q*m)+vec2(0.1, 4.741)), fbm((p+q*m)+vec2(1.952, 7.845))); 
            m /= mm;
            return fbm(p+r*m);
        }

        float fbm( vec2 x) {
            perlin_noise noise(1998);
            float h = 0.0;

            for (float i=1.0;i<10.0;i++) {
                h+=noise.eval(x*pow(1.6, i))*0.9*pow(0.6, i);
            }
            return h;
        }

        


};

class liquid {
    public:
        float eval(vec2 p, float frequency) {
            float val = pow(fbm(warp(p * frequency, 8.0, 4) + 20.0, 5), 2.2);
            return val * 2;
        }

        float fbm(vec2 co, int layers) {
            float val = 0.0;
            float factor = 1.0;
            for (int n = 0; n < layers; n++) {
                val += (noise.eval(co * factor) - 0.5)/factor;
                factor *= 2.0;
            }
            return (val + 1.0) * 0.5;
        }

        vec2 warp(vec2 co, float intensity, int detail) {
            return co + vec2((fbm(co, detail) - .5) * intensity, (fbm(co + vec2(12342.145, -2340.769), detail) - .5) * intensity);
        }
    private:
        perlin_noise noise = perlin_noise(1998);
        
};

class turbunlent_liquid {
    public:
        float eval(vec2 p, float frequency) {
            float val = pow(fbm(warp(p * frequency, 8.0, 4) + 20.0, 5), 2.2) + pow(fbm2(warp(p * frequency, 8.0, 4) + 20.0, 5), 2.2) ;
            return val * 2;
        }

        float fbm(vec2 co, int layers) {
            float val = 0.0;
            float factor = 1.0;
            for (int n = 0; n < layers; n++) {
                val += (noise.eval(co * factor,1,1) - 0.5)/factor;
                factor *= 2.0;
            }
            return (val + 1.0) * 0.5;
        }

        float fbm2(vec2 co, int layers) {
            float val = 0.0;
            float amp = 1.0;
            float frequency = 0.01;
            float lacunarity = 1.8;
            vec2 p = co * frequency;
            for (int n = 0; layers < layers; n++) {
                val += std::fabs(2 *noise.eval(p, 1,1) - 1) * amp;
                amp *= 2.0;
                p  *= lacunarity;
            }
            return val * 0.2;
        }

        vec2 warp(vec2 co, float intensity, int detail) {
            return co + vec2((fbm(co, detail) - .5) * intensity, (fbm(co + vec2(12342.145, -2340.769), detail) - .5) * intensity);
        }

    private:
        voronoise noise;
        
};




#endif