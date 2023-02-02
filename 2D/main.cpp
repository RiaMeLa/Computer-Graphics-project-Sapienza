
#include "vec2.h"
#include "image.h"
#include "noiselib.h"
#include <cstdio> 
#include <random> 
#include <functional> 
#include <iostream> 
#include <fstream> 
#include <cmath> 





int main() {
    int w = 512;
    int h = 512;


/*##############VORONOISE###################

    voronoise vor;
    vec3* img_buf = new vec3[h * w];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p = p / 512;
            p = p*(16) + 14;
            float d = vor.eval(p, 1, 0.75);
            img_buf[i * w + j] = vec3(d, d, d);
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*##########################################*/

/*###############SMOOTH VORONOI############
    smooth_voronoi vor;
    vec3* img_buf = new vec3[h * w];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p = p / 512;
            p = p*(4);
            float d = vor.eval(p);
            img_buf[i * w + j] = vec3(d, d, d);
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*#########################################*/

/*###########COLORED VORONOI############
    voronoi_colored vor2;
    vec3* img_buf = new vec3[w*h];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p = p / 512;
            p = (p*6) + 14.0;
            img_buf[i * w + j] = vor2.eval(p);
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*###########COLORED VORONOI############*/

/*##################PERLIN NOISE########################
    vec3* img_buf = new vec3[w * h];
    float noise_map [w * h];
    perlin_noise noise(1998);

    float frequency  = 0.02;
    float lacunarity = 1.8;
    float gain       = 0.35;
    int   octaves    = 1;
    float max_noise_val = 0.0;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p *= frequency;
            float amplitude = 1;
            float val = 0.0;
            noise_map[i * w + j] = 0.0;
            for (int l = 0; l < octaves; ++l) {
                //fractal sum
                noise_map[i * w + j] += noise.eval(p) * amplitude;
                p *= lacunarity;
                amplitude *= gain;
            }
            //if (noise_map[i * w + j] > max_noise_val) max_noise_val = noise_map[i * w + j];
        }
    }
    //for (int i = 0; i < w * h; ++i) noise_map[i] /= max_noise_val;
    img_buf = noisemap_to_grayscale(noise_map, w, h);
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*###################VALUE NOISE##########################*/


/*################## MARBLE ########################
    vec3* img_buf = new vec3[w * h];
    perlin_noise noise(1998);

    float frequency  = 0.02;
    float lacunarity = 1.8;
    float gain       = 0.35;
    int octaves      = 5;

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            //fbm
            p *= frequency;
            float amplitude = 1.0;
            float noise_val = 0.0;
            for (int l = 0; l < octaves; ++l) {
                noise_val += noise.eval(p) * amplitude + (std::fabs(2 *noise.eval(p) - 1) * amplitude/1.5);
                p *= lacunarity;
                amplitude *= gain;
            }
            img_buf[i * w + j] = vec3(  (sin((j + noise_val * 100) * 2 * M_PI / 200.f) + 1) / 2.f,
                                        (sin((j + noise_val * 100) * 2 * M_PI / 200.f) + 1) / 2.f,
                                        (sin((j + noise_val * 100) * 2 * M_PI / 200.f) + 1) / 2.f );     
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*#######################################################*/

/*################## ABSTRACT 1 ########################
    vec3* img_buf = new vec3[w * h];
    perlin_noise noise(8);
    voronoise vor;
    float frequency  = 0.02;
    float lacunarity = 1.8;
    float gain       = 0.35;
    int octaves      = 5;
    float max_noise_val = 0.0;

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p = vec2(i, j) * frequency;
            //fbm
            float amplitude = 1.0;
            float noise_val = 0.0;
            for (int l = 0; l < octaves; ++l) {
                noise_val += vor.eval(vec2(p.x + noise.eval(p), p.y + noise.eval(p)), 1, 1) * amplitude;
                p *= lacunarity;
                amplitude *= gain;
            }

            float noise_val2 = 0.0;   
            p = vec2(i, j) * 0.01;
            amplitude = 1;
            for (int l = 0; l < octaves; ++l) {
                noise_val2 += std::fabs(2 *noise.eval(vec2(p.x + vor.eval(p, 1, 1), p.y + vor.eval(p, 1, 1))) - 1) * amplitude;
                p *= 2.5;
                amplitude *= 0.15;
            }
            if(max_noise_val < noise_val2) max_noise_val = noise_val2;
            img_buf[i * w + j] = vec3(noise_val, noise_val2, 0);
        }
    }
    for (int i = 0; i < h; ++i) 
    for (int j = 0; j < w; ++j) {
        vec3 p = img_buf[i * w + j];
        p = vec3( (sin((1 + (p.x + p.y/max_noise_val)* 155) * 2 * M_PI / 200.f) + 1) / 2.f,
                  (sin((1 + (p.x + p.y/max_noise_val)* 200) * 2 * M_PI / 200.f) + 1) / 2.f,
                  (sin((1 + (p.x + p.y/max_noise_val)* 200) * 2 * M_PI / 200.f) + 1) / 2.f);
        img_buf[i * w + j] = p;
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*#######################################################*/

/*#######################VORONOI#########################
    vec3* img_buf = new vec3[w * h];
    smooth_voronoi vor;
    vec3* img_buf = new vec3[h * w];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p = p / std::max(w, h);
            p = p*(3);
            float d = vor.eval(p);
            img_buf[i * w + j] = vec3(d, d, d);
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*#######################################################*/

/*###################VORONOI DISTANCES##################
    vec3* img_buf = new vec3[w * h];
    voronoi_distances vor;
    
    for (int i = 0; i < h; ++i)
        for(int j = 0; j < w; ++j) {
            vec2 p(i, j);
            p = p / std::max(w, h);
            vec3 c = vor.eval(p * 8.0);

            vec3 col = c.x*(0.5 + 0.5*sin(64.0*c.x))*vec3(1.0, 1.0, 1.0);

            col = lerp(vec3(1.0,0.6,0.0), col, smoothstep( 0.04, 0.07, c.x ));

            float dd = sqrt(dot(vec2(c.y, c.z), vec2(c.y, c.z)));
            //col = lerp( vec3(1.0,0.6,0.01), col, smoothstep( 0.0, 0.12, dd));
            img_buf[i * w  + j] = col;
        }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*######################################################*/

/*###################TURBULENCE VARIATION#############################
    vec3* img_buf = new vec3[w * h];
    float noise_map[h * w];
    
    voronoise noise;
    voronoi vor;
    float frequency  = 0.007;
    float lacunarity = 3.8;
    float gain       = 0.35;
    int   octaves    = 5;
    float max_noise_val = 0.0;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p = vec2(i,j) * frequency;
            float amplitude = 1;
            noise_map[i * w + j] = noise.eval(vec2(p.x + vor.eval(p), p.y + vor.eval(p)), 0, 1);
            for (int l = 0; l < octaves; ++l) {
                //fractal sum
                noise_map[i * w + j] += std::fabs(2 *noise.eval(vec2(p.x + vor.eval(p), p.y + vor.eval(p)), 0, 1) - 1) * amplitude;
                p *= lacunarity;
                amplitude *= gain;
            }
            if (noise_map[i * w + j] > max_noise_val) max_noise_val = noise_map[i * w + j];
        }
    }
    for (int i = 0; i < w * h; ++i) noise_map[i] /= max_noise_val;
    
    
    img_buf = noisemap_to_grayscale(noise_map, w, h);
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*##############################################################*/

/*###################TURBULENCE#############################
    vec3* img_buf = new vec3[w * h];
    float noise_map[h * w];
    
    perlin_noise noise(1009);
    float frequency  = 0.02;
    float lacunarity = 1.8;
    float gain       = 0.35;
    int   octaves    = 5;
    float max_noise_val = 0.0;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            vec2 p = vec2(i,j) * frequency;
            float amplitude = 1;
            noise_map[i * w + j] = 0.0;
            for (int l = 0; l < octaves; ++l) {
                //fractal sum
                noise_map[i * w + j] += std::fabs(2 *noise.eval(p) - 1) * amplitude;
                p *= lacunarity;
                amplitude *= gain;
            }
            if (noise_map[i * w + j] > max_noise_val) max_noise_val = noise_map[i * w + j];
        }
    }
    for (int i = 0; i < w * h; ++i) noise_map[i] /= max_noise_val;
    
    
    img_buf = noisemap_to_grayscale(noise_map, w, h);
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*##############################################################*/

/*###############LIQUID###############*/
    vec3* img_buf = new vec3[w * h];
    float detail = 4; 
    float intensity = 8;
    float frequency = 0.01;
    turbunlent_liquid noise;
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            float val = noise.eval(vec2(i,j), frequency);
            img_buf[i * w + j] = vec3(val, val, val);
        }
    }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*###############################*/

/*##########CLOUDS###########
    vec3* img_buf = new vec3[w * h];
    cloudy noise;
    for (int i = 0; i < h; ++i) 
        for (int j = 0; j < w; ++j) {
            vec2 p = vec2(i,j);
            float col = noise.warp(p*0.004, 12.0+noise.fbm(p*0.005)*16.0);
            float y = pow(1.0-p.y/w, 2.0);
            vec3 color = lerp(vec3(0.2+0.3*y, 0.4+0.2*y, 1.0), vec3(1,1,1), smoothstep(0.5, 1.0, col));
            img_buf[i * w + j] = color;
        }
    image img(w, h, img_buf);
    print_image(std::cout, img);
    delete[] img_buf;
/*#############################*/




}