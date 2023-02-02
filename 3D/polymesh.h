#include "vec3.h"

class polymesh {
    public:
        polymesh() : vertices(nullptr), st(nullptr), normals(nullptr) {}
        ~polymesh() {
            if (vertices) delete[] vertices;
            if (st) delete[] st;
            if (normals) delete[] normals;
        }
    public:
        vec3* vertices;
        vec3* st;
        vec3* normals;
        int* face_array;
        int* vertices_array;
        int  num_vertices;
        int  num_faces;
};

polymesh* create_polymesh(
    int width = 1,
    int height = 1,
    int subdiv_width = 1,
    int subdiv_height = 1)
{
    polymesh* poly = new polymesh;
    poly->num_vertices = (subdiv_width + 1) * (subdiv_height + 1);
    poly->vertices = new vec3[poly->num_vertices];
    poly->normals = new vec3[poly->num_vertices];
    poly->st = new vec3[poly->num_vertices];
    float inv_subdiv_w = 1.f / subdiv_width;
    float inv_subdiv_h = 1.f / subdiv_height;
    vec3 v;
    for (int i = 0; i <= subdiv_height; ++i)
        for (int j = 0; j <= subdiv_width; ++j) {
            v = vec3(1-i, 1-j, -1);;
            //v = vec3(width * (j * inv_subdiv_w - 0.5), height * (i * inv_subdiv_h - 0.5), 0);
            poly->vertices[i * (subdiv_height + 1) + j] = v;
            poly->st[i * (subdiv_width + 1) + j] = vec3(j * subdiv_width, i * subdiv_height, 0); //being vec2
        }

    poly->num_faces = subdiv_width * subdiv_height;
    poly->face_array = new int[poly->num_faces];
    for (int i = 0; i < poly->num_faces; ++i)
        poly->face_array[i] = 4;

    poly->vertices_array = new int[4 * poly->num_faces];
    for (int i = 0, k = 0; i < subdiv_height; ++i)
        for (int j = 0; j < subdiv_width; ++j) {
            poly->vertices_array[k] = i * (subdiv_width + 1) + j;
            poly->vertices_array[k + 1] = i * (subdiv_width + 1) + j + i;
            poly->vertices_array[k + 2] = (i + 1) * (subdiv_width + 1) + j + 1;
            poly->vertices_array[k + 3] = (i + 1) * (subdiv_width + 1) + j;
            k += 4;
        }
    
    return poly;
}
