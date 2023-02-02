#ifndef HITTABLE_H
#define HITTABLE_H
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
#include "ray.h"
#include <math.h>
class material;

const float epsilon = 1e-8;

struct hit_record {
    point3 p;
    vec3 normal;
    shared_ptr<material> mat_ptr;
    float t;
    float u;
    float v;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal :-outward_normal;
    }
};


class hittable {
    public:
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};

class triangle : public hittable {
    public:
        triangle(point3 point0, point3 point1, point3 point2, shared_ptr<material> m) : p0(point0), p1(point1), p2(point2), mat_ptr(m) {}
         
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override {
            //compute plane normals
            vec3 p0p1 = p1 - p0;
            vec3 p0p2 = p2 - p0;
            vec3 n = unit_vector(cross(p0p1, p0p2));
            //float area_squared = dot(n, n);

            /*##finding ray-plane point of intersection##*/

            float n_dot_raydir = dot(n, r.direction());
            if (fabs(n_dot_raydir) < epsilon) return false; //if parallel there's no intersection
            //compute D param of ax + by + cz + d = 0 equation, 1 point is 
            //enough since they're complanar
            float d = -dot(n, p0);
            //compute t o fray equation r(t) = o+dir(t)
            rec.t = -(dot(n, r.origin()) + d) / n_dot_raydir;
            if (rec.t < t_min || t_max < rec.t) return false;
            //check if behind ray
            //if (rec.t < 0) return false;
            //compute intersection point
            vec3 p = r.origin() + rec.t * r.direction();
            
            /*##inside-outside  test###*/
            //vector perpendicular to triangles plane
            vec3 c;
            //edge0
            vec3 edge0 = p1 - p0;
            vec3 vp0 = p - p0;
            c = cross(edge0, vp0);
            if (dot(n, c) < 0) return false; // p is on the right
        
            //edge1
            vec3 edge1 = p2 - p1;
            vec3 vp1 = p - p1;
            c = cross(edge1, vp1);
            if (dot(n, c) < 0) return false;

            //edge2
            vec3 edge2 = p0 - p2;
            vec3 vp2 = p - p2;
            c = cross(edge2, vp2);
            if (dot(n, c) < 0) return false;

            rec.p = p;
            rec.mat_ptr = mat_ptr;
            rec.set_face_normal(r, n);
            return true; // on the triangle if on the left of each edge
        }   
    private:
        point3 p0, p1, p2;
        shared_ptr<material> mat_ptr;
};

class rectangle : public hittable {
    public:
        rectangle(point3 point0, point3 point1, point3 point2, point3 point3, shared_ptr<material> m) : p0(point0), p1(point1), p2(point2), mat_ptr(m) {}
         
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override {
            //compute plane normals
            vec3 p0p1 = p1 - p0;
            vec3 p0p2 = p2 - p0;
            vec3 n = unit_vector(cross(p0p1, p0p2));
            //float area_squared = dot(n, n);

            /*##finding ray-plane point of intersection##*/

            float n_dot_raydir = dot(n, r.direction());
            if (fabs(n_dot_raydir) < epsilon) return false; //if parallel there's no intersection
            //compute D param of ax + by + cz + d = 0 equation, 1 point is 
            //enough since they're complanar
            float d = -dot(n, p0);
            //compute t o fray equation r(t) = o+dir(t)
            rec.t = -(dot(n, r.origin()) + d) / n_dot_raydir;
            if (rec.t < t_min || t_max < rec.t) return false;
            //check if behind ray
            //if (rec.t < 0) return false;
            //compute intersection point
            vec3 p = r.origin() + rec.t * r.direction();
            
            /*##inside-outside  test###*/
            //vector perpendicular to triangles plane
            vec3 c;
            //edge0
            vec3 edge0 = p1 - p0;
            vec3 vp0 = p - p0;
            c = cross(edge0, vp0);
            if (dot(n, c) < 0) return false; // p is on the right
        
            //edge1
            vec3 edge1 = p2 - p1;
            vec3 vp1 = p - p1;
            c = cross(edge1, vp1);
            if (dot(n, c) < 0) return false;

            //edge2
            vec3 edge2 = p3 - p2;
            vec3 vp2 = p - p2;
            c = cross(edge2, vp2);
            if (dot(n, c) < 0) return false;

            //edge2
            vec3 edge3 = p0 - p3;
            vec3 vp3 = p - p3;
            c = cross(edge3, vp3);
            if (dot(n, c) < 0) return false;

            rec.p = p;
            rec.mat_ptr = mat_ptr;
            rec.set_face_normal(r, n);
            return true; // on the triangle if on the left of each edge
        }   
    private:
        point3 p0, p1, p2, p3;
        shared_ptr<material> mat_ptr;
};


#endif
