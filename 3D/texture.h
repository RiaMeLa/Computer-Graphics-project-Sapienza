#ifndef TEXTURE_H
#define TEXTURE_H
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

#include "noiselib.h"
#include "vec3.h"
#include <iostream>


class texture  {
    public:
        virtual color value(float u, float v, const vec3& p) const = 0;
};


class solid_color : public texture {
    public:
        solid_color() {}
        solid_color(color c) : color_value(c) {}

        solid_color(float red, float green, float blue)
          : solid_color(color(red,green,blue)) {}

        virtual color value(float u, float v, const vec3& p) const override {
            return color_value;
        }

    private:
        color color_value;
};


class checker_texture : public texture {
    public:
        checker_texture() {}

        checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
            : even(_even), odd(_odd) {}

        checker_texture(color c1, color c2)
            : even(make_shared<solid_color>(c1)) , odd(make_shared<solid_color>(c2)) {}

        virtual color value(float u, float v, const vec3& p) const override {
            auto sines = sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
            if (sines < 0)
                return odd->value(u, v, p);
            else
                return even->value(u, v, p);
        }

    public:
        shared_ptr<texture> odd;
        shared_ptr<texture> even;
};


class perlin_noise_texture : public texture {
    public:
        perlin_noise_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(noise.eval(p * 8));
        }

    public:
        perlin_noise noise;
};

class value_noise_texture : public texture {
    public:
        value_noise_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(noise.eval(p * 8));
        }

    public:
        value_noise noise;
};

class voronoi_texture : public texture {
    public:
        voronoi_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(vor.eval(p * 8));
        }

    public:
        voronoi vor;
};

class voronoi_colored_texture : public texture {
    public:
        voronoi_colored_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*2*(vor.eval(p * 0.3));
        }

    public:
        voronoi_colored vor;
};

class smooth_voronoi_texture : public texture {
    public:
        smooth_voronoi_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(vor.eval(p * 8));
        }

    public:
        smooth_voronoi vor;
};

class voronoise_texture : public texture {
    public:
        voronoise_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(vor.eval(p * 10, 1, 1));
        }

    public:
        voronoise vor;
};

class voronoi_distances_texture : public texture {
    public:
        voronoi_distances_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            return color(1,1,1)*(vor.eval(p * 2));
        }

    public:
        voronoi_distances vor;
};

class marble_texture : public texture {
    public:
        marble_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            float frequency  = 1.0;
            float lacunarity = 2.0;
            float gain       = 0.5;
            int octaves      = 5;
            float max_noise_val = 0.0;

            vec3 ps = p * frequency * 5;
            float amplitude = 1.0;
            float noise_val = 0.0;

            for (int l = 0; l < octaves; ++l) {
                noise_val += noise.eval(ps) * amplitude + std::fabs(1.5 * noise.eval(ps) - 1) * amplitude;
                ps *= lacunarity;
                amplitude *= gain;
            }



            
            return color(1, 1, 1)* 0.5 * (1 + sin(5 * p.z() + 10*noise_val));
            //return color(1, 1, 1) * (sin((p.z() + noise_val * 100) * 2 * M_PI / 200.f) + 1) / 2.f;
        }

    public:
        perlin_noise noise;
};

class abstract1_texture : public texture {
    public:
        abstract1_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            float frequency  = 0.01;
            float lacunarity = 1.8;
            float gain       = 0.35;
            int octaves      = 5;
            float max_noise_val = 0.0;

            vec3 ps = p * frequency;
            float noise_val = 0.0;
            float amplitude = 1.0;
            for (int l = 0; l < octaves; ++l) {
                noise_val += noise.eval(vec3(p.x() + noise.eval(p), p.y() + noise.eval(p), p.z() + noise.eval(p))) * amplitude;
                ps *= lacunarity;
                amplitude *= gain;
            }

            float noise_val2 = 0.0;   
            ps = p * 0.008;
            amplitude = 1;
            for (int l = 0; l < octaves; ++l) {
                noise_val2 += std::fabs(2 *noise.eval(vec3(p.x() + noise.eval(p), p.y() + noise.eval(p), p.z() + noise.eval(p))) - 1) * amplitude;
                ps *= 2.5;
                amplitude *= 0.15;
            }           
            ps = vec3 ( (sin((1 + (noise_val + noise_val2)* 170) * 2 * M_PI / 200.f) + 1) / 2.f,
                        (sin((1 + (noise_val + noise_val2)* 200) * 2 * M_PI / 200.f) + 1) / 2.f,
                        (sin((1 + (noise_val + noise_val2)* 200) * 2 * M_PI / 200.f) + 1) / 2.f );

            return color(1, 1, 1) * ps;
        }

    public:
        perlin_noise noise;
};

class turbulence_texture : public texture {
    public:
        turbulence_texture() {}

        virtual color value(float u, float v, const vec3& p) const override {
            // return color(1,1,1)*0.5*(1 + noise.turb(scale * p));
            // return color(1,1,1)*noise.turb(scale * p);
            float frequency  = 0.9;
            float lacunarity = 1.8;
            float gain       = 0.35;
            int octaves      = 5;
            float max_noise_val = 0.0;

            vec3 ps = p * frequency;
            float amplitude = 1.0;
            float noise_val = 0.0;

            for (int l = 0; l < octaves; ++l) {
                noise_val += (2 * noise.eval(ps) - 1) * amplitude;
                ps *= lacunarity;
                amplitude *= gain;
            }
            return color(1, 1, 1)* fabs(noise_val);
        }

    public:
        perlin_noise noise;
};


#endif
