#ifndef PERLIN_H
#define PERLIN_H
//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"
#include "utils.h"
#include "remappers.h"
#include "utils.h"

class value_noise {
    public:
        value_noise() {
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

        float eval(const vec3& p) const {
            int xi = std::floor(p.x());
            int yi = std::floor(p.y());
            int zi = std::floor(p.z());

            float tx = p.x() - xi;
            float ty = p.y() - yi;
            float tz = p.z() - zi;

            float sx = smoothstep(tx);
            float sy = smoothstep(ty);
            float sz = smoothstep(tz);

            int xi0 = xi & max_size_mask;
            int xi1 = (xi0 + 1) & max_size_mask;
            int yi0 = yi & max_size_mask;
            int yi1 = (yi0 + 1) & max_size_mask;
            int zi0 = zi & max_size_mask;
            int zi1 = (zi0 + 1) & max_size_mask;

            float c000 = r[permutation_table[permutation_table[permutation_table[xi0] + yi0] + zi0]];
            float c100 = r[permutation_table[permutation_table[permutation_table[xi1] + yi0] + zi0]];            
            float c010 = r[permutation_table[permutation_table[permutation_table[xi0] + yi1] + zi0]];
            float c110 = r[permutation_table[permutation_table[permutation_table[xi1] + yi1] + zi0]];
            float c001 = r[permutation_table[permutation_table[permutation_table[xi0] + yi0] + zi1]];
            float c101 = r[permutation_table[permutation_table[permutation_table[xi1] + yi0] + zi1]];
            float c011 = r[permutation_table[permutation_table[permutation_table[xi0] + yi1] + zi1]];
            float c111 = r[permutation_table[permutation_table[permutation_table[xi1] + yi1] + zi1]];
    
            // linearly interpolate values along the x axis
            float a = lerp(c000, c100, sx); 
            float b = lerp(c010, c110, sx);
            float c = lerp(c001, c101, sx);
            float d = lerp(c011, c111, sx);

            float e = lerp(a, b, sy);
            float f = lerp(c, d, sy);

            return lerp(e, f, sz);
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
        perlin_noise () {
            for (int i = 0; i < table_size; ++i) { 
                float theta = acos(2 * random_float(0, 1) - 1);
                float phi = 2 * random_float(0, 1) * M_PI;

                float x = cosf(phi) * sinf(theta);
                float y = sinf(phi) * sinf(theta);
                float z = cosf(theta);
                gradients[i] = vec3(x, y, z);
                permutation_table[i] = i; 
            }
            //shuffle permutation table
            for (int i = 0; i < table_size; ++i) 
                std::swap(permutation_table[i], permutation_table[rand() & mask_size]); 
            // extend the permutation table in the index range [256:512]
            for (int i = 0; i < table_size; ++i) 
                permutation_table[table_size + i] = permutation_table[i]; 
        }

        float eval(const vec3& p) const {  
            /*floor returns lowest int if input is signed*/
            int xi  = std::floor(p.x());
            int yi  = std::floor(p.y());
            int zi  = std::floor(p.z());

            int xi0 = xi & mask_size;
            int yi0 = yi & mask_size;
            int zi0 = zi & mask_size;
            int xi1 = (xi0 + 1) & mask_size;
            int yi1 = (yi0 + 1) & mask_size; 
            int zi1 = (zi0 + 1) & mask_size; 

            float tx = p.x() - xi;
            float ty = p.y() - yi;
            float tz = p.z() - zi;

            float sx = quintic(tx);
            float sy = quintic(ty);
            float sz = quintic(tz);

            //gradients instead of floats at corners
            vec3 c000 = gradients[permutation_table[permutation_table[permutation_table[xi0] + yi0] + zi0]];
            vec3 c100 = gradients[permutation_table[permutation_table[permutation_table[xi1] + yi0] + zi0]];            
            vec3 c010 = gradients[permutation_table[permutation_table[permutation_table[xi0] + yi1] + zi0]];
            vec3 c110 = gradients[permutation_table[permutation_table[permutation_table[xi1] + yi1] + zi0]];
            vec3 c001 = gradients[permutation_table[permutation_table[permutation_table[xi0] + yi0] + zi1]];
            vec3 c101 = gradients[permutation_table[permutation_table[permutation_table[xi1] + yi0] + zi1]];
            vec3 c011 = gradients[permutation_table[permutation_table[permutation_table[xi0] + yi1] + zi1]];
            vec3 c111 = gradients[permutation_table[permutation_table[permutation_table[xi1] + yi1] + zi1]];
            

            //vectors from grid points to p
            float x0 = tx; 
            float y0 = ty; 
            float z0 = tz; 
            float x1 = tx - 1;
            float y1 = ty - 1;
            float z1 = tz - 1;

            vec3 p000 = vec3(x0, y0, z0); 
            vec3 p100 = vec3(x1, y0, z0); 
            vec3 p010 = vec3(x0, y1, z0); 
            vec3 p110 = vec3(x1, y1, z0); 
    
            vec3 p001 = vec3(x0, y0, z1); 
            vec3 p101 = vec3(x1, y0, z1); 
            vec3 p011 = vec3(x0, y1, z1); 
            vec3 p111 = vec3(x1, y1, z1);

            //lerp
            float a = lerp(dot(c000, p000), dot(c100, p100), sx); 
            float b = lerp(dot(c010, p010), dot(c110, p110), sx); 
            float c = lerp(dot(c001, p001), dot(c101, p101), sx); 
            float d = lerp(dot(c011, p011), dot(c111, p111), sx);

            float e = lerp(a, b, sy); 
            float f = lerp(c, d, sy); 
    
            return (lerp(e, f, sz) + 1) * 0.5;  
        }
    private:
        int seed;
        static const int table_size = 256; 
        static const int mask_size = table_size - 1; 
        vec3 gradients[table_size]; 
        unsigned permutation_table[table_size * 2]; 
};

class voronoi {
    public:
        voronoi() {}

        float eval(const vec3& p) const {
            vec3 n = vec3(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
            vec3 f = fract(p);

            float res = 80.0;
            for (int x = -1; x <= 1; ++x)
                for (int y = -1; y <= 1; ++y)
                    for (int z = -1; z <= 1; ++z) {
                        vec3 b = vec3(float(x), float(y), float(z));
                        //vec2 k = hash(b + n);
                        vec3 k = n33(n + b);
                        vec3 r = b - f + k;
                        //r = r + k;
                        float d = dot(r, r);
                        res = std::min(res, d);
                }
            return sqrt(res);
        }

    private:
        vec3 n33(const vec3& p) const {
            vec3 a = fract(p * vec3(564.391,309.444,666.777));
            a +=  dot(a, a*5.63);
            return fract(vec3(a.x()*a.y(), a.y()*a.z(), a.x()*a.z()));
        }
};

class voronoi_colored {
    public:
        voronoi_colored() {}

        vec3 eval(const vec3& p) const {
            vec3 n = vec3(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
            vec3 f = fract(p);

            float res = 80.0;
            vec3 col;
            for (int x = -1; x <= 1; ++x)
                for (int y = -1; y <= 1; ++y)
                    for (int z = -1; z <= 1; ++z) {
                        vec3 b = vec3(float(x), float(y), float(z));
                        //vec2 k = hash(b + n);
                        vec3 k = n33(n + b);
                        vec3 r = b - f + k;
                        //r = r + k;
                        float d = dot(r, r);
                        if (d < res) {
                            res = d;                         //k should be b + n
                            vec3 v = 0.5 + 0.5 * (hash1(dot(b + n,vec3(7.0,13.0, 3.18)))*1.5 + 1.8 + vec3(k.x() + 2.7, k.x() +  1.9, k.x() + 1.7));
                            col = vec3(sinf(v.x()) ,sinf(v.y()) ,sinf(v.z()));
                            col = col * col;
                        }
                }
            return col;
        }

    private:
        vec3 n33(const vec3& p) const {
            vec3 a = fract(p * vec3(564.391,309.444,666.777));
            a +=  dot(a, a*5.63);
            return fract(vec3(a.x()*a.y(), a.y()*a.z(), a.x()*a.z()));
        }
        float hash1(const float n ) const { return fract(sin(n)*43758.5453); }
};

class smooth_voronoi {
    public:
        smooth_voronoi() {}

        float eval(const vec3& p) const {
            vec3 n = vec3(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
            vec3 f = fract(p);

            float res = 0.0;
            for (int x = -1; x <= 1; ++x)
                for (int y = -1; y <= 1; ++y)
                    for (int z = -1; z <= 1; ++z) {
                        vec3 b = vec3(float(x), float(y), float(z));
                        //vec2 k = hash(b + n);
                        vec3 k = n33(n + b);
                        vec3 r = b - f + k;
                        //r = r + k;
                        float d = dot(r, r);
                        res += 1.0/pow( d, 12.0 );
                    }
            return pow( 1.0/res, 1.0/16.0 );
        }

    private:
        vec3 n33(const vec3& p) const {
            vec3 a = fract(p * vec3(564.391,309.444,666.777));
            a +=  dot(a, a*5.63);
            return fract(vec3(a.x()*a.y(), a.y()*a.z(), a.x()*a.z()));
        }
};

class voronoise {
    public: 
        voronoise() {}

        float eval(const vec3& p, float u, float v) const {
            vec3 n = vec3(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
            vec3 f = fract(p);

            float sharpness = 1.0 + 63.0 * pow(1.0-v, 4.0);

            float value = 0.0;
            float accum = 0.0;
            for (int x = -2; x <= 2; ++x)
                for (int y = -2; y <= 2; ++y)
                    for(int z = -2; z <= 2; ++z) {
                    vec3 curr = vec3(float(x), float(y), float(z));
                    vec3 cen = n33(curr + n, u);
                    vec3 r = curr - f + cen;
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
        float rand(const vec3 &p) const {
            return fract(sin(dot(p,vec3(419.2,371.9, 287.2))) * 833458.57832);
        }

        vec3 n33(const vec3& p, float u) const {
            vec3 a = fract(p * vec3(564.391,309.444,666.777));
            a +=  dot(a, a*5.63);
            return fract(vec3(a.x()*a.y(), a.y()*a.z(), a.x()*a.z()) * u);
        }
};

class voronoi_distances {
    public:
        vec3 eval(const vec3& p) const {
            /*As first we search which cell contains the closest center 
            to our point p*/
            vec3 n = vec3(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
            vec3 f = fract(p);

            vec3 mg, mr;
            float md = 8.0;
            for(int x=-1; x<=1; x++)
            for(int y=-1; y<=1; y++)
            for(int z=-1; z<=1; z++) {
                vec3 curr = vec3(float(x), float(y), float(z));
                vec3 cen = n33(curr + n);
                vec3 r = curr + cen - f;
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
            for(int x=-2; x<=2; x++)
            for(int y=-2; y<=2; y++)
            for(int z=-2; z<=2; z++) {
                vec3 g = mg + vec3(float(x), float(y), float(z));
                vec3 cen = n33(n + g);
                vec3 r = g + cen - f;
                if( dot(mr-r,mr-r)>0.00001 )
                    md = std::min(md, dot(0.5 * (mr + r), unit_vector(r - mr)));
            
            }
            vec3 col = md*(0.5 + 0.5*sin(64.0*md))*vec3(1.0, 1.0, 1.0);
            col = lerp(vec3(1.0,0.6,0.0), col, smoothstep( 0.04, 0.07, md));
            return col;
        }

        vec3 n33(const vec3& p) const {
            vec3 a = fract(p * vec3(564.391,309.444,666.777));
            a +=  dot(a, a*5.63);
            return fract(vec3(a.x()*a.y(), a.y()*a.z(), a.x()*a.z()));
        }
};

#endif
