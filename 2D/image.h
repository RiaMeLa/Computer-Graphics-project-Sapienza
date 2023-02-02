#ifndef IMAGE_H
#define IMAGE_H

#include "vec3.h"

#include <math.h>
#include "utils.h"
#include "image.h"
#include <cstdio> 
#include <random> 
#include <functional> 
#include <iostream> 
#include <fstream> 
#include <cmath> 

class image {
    public: 
        image(int w, int h) : width(w), height(h) {
            img_buffer = (vec3*) malloc(h*w*sizeof(vec3));
        }
        image(int w, int h, vec3* img_buf) : width(w), height(h), img_buffer(img_buf) {}

        int h() {
            return height;
        }

        int w() {
            return width;
        }
        
    public:
        int height;
        int width;
        vec3* img_buffer;
};

void print_image(std::ostream &out, image &img) {
    vec3* buf = img.img_buffer;
    out << "P3\n" << img.w() << ' ' << img.h() << "\n255\n";
    for (int i = 0; i < img.h(); ++i) {
        for (int j = 0; j < img.w(); ++j) {
            float r = img.img_buffer[i * img.w() + j].x;
            float g = img.img_buffer[i * img.w() + j].y;
            float b = img.img_buffer[i * img.w() + j].z;
            out << static_cast<int>(255 * r) << ' '
                << static_cast<int>(255 * g) << ' '
                << static_cast<int>(255 * b) << '\n';
        }
    }

}

vec3* noisemap_to_grayscale(float* noise_map, int w, int h) {
    vec3* gray_img = new vec3[w * h];
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
                float r = noise_map[i * w + j];
                float g = noise_map[i * w + j];
                float b = noise_map[i * w + j];
                gray_img[i * w + j] = vec3(r, g, b);
        }
    }
    return gray_img;
}

#endif