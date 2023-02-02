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

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "texture.h"
#include "polymesh.h"
#include <iostream>


color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}


int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 10;

    // World

    hittable_list world;
    auto perlintext = perlin_noise_texture();
    auto voronoicolored = voronoi_colored_texture();
    auto vorotext = voronoi_texture();
    auto smoothvortext = smooth_voronoi_texture();
    auto vorodist = voronoi_distances_texture();
    auto voronoisytext = voronoise_texture();
    auto marbletext = marble_texture();
    auto turbulenttext = turbulence_texture();
    auto abstract = abstract1_texture();
    auto solidcol = solid_color(1,0,0);

    polymesh* poly = create_polymesh();
//  USE THESE TO CHANGE SPHERE TEXTURE
    auto red = make_shared<lambertian>(&solidcol);
    auto perlin = make_shared<lambertian>(&perlintext);
    auto voronoise = make_shared<lambertian>(&voronoisytext);
    auto abstractive = make_shared<lambertian>(&abstract);
    auto marble = make_shared<lambertian>(&marbletext);
    auto turbulent = make_shared<lambertian>(&turbulenttext);
    auto voronoidist = make_shared<lambertian>(&vorodist);
    auto vorocolored = make_shared<lambertian>(&voronoicolored);
    auto voronoi = make_shared<lambertian>(&vorotext);
    auto smoothvor = make_shared<lambertian>(&smoothvortext);

    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, red));
    //REPLACE NAME OF TEXTURE--------------------------------------------------------->HERE
    world.add(make_shared<triangle>(point3(0,0,-1), point3(5,0,-1), point3(2.5,3,-1), perlin));
    world.add(make_shared<sphere>(point3(1,2,-1), 2, voronoi));
    world.add(make_shared<sphere>(point3(6,1,4), 1, marble));
    world.add(make_shared<sphere>(point3(6,1,-1.5), 1, turbulent));
    world.add(make_shared<sphere>(point3(-4,2.5,2.5), 2.5, voronoise));


    // Camera

    point3 lookfrom(13,2,3);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;

    camera cam(lookfrom, lookat, vup, 50, aspect_ratio, aperture, dist_to_focus);

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0,0,0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_float()) / (image_width-1);
                auto v = (j + random_float()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}


