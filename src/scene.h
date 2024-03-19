#pragma once

#include "lib/image.h"
#include "lib/ray_tracer.h"
#include "loader.h"
#include "materials/diffuse.h"

struct camera {
    Vec3f position;
    Vec3f look_dir;
    Vec3f up;
    float fov;
    float aspect;

    camera(){
        position = {0.25, -0.25, 0.5};
        look_dir = normalize(Vec3f(0.0, 0.0, 0) - position);
        up = {1, 0, 0};
        fov = 80;
        aspect = 1;
    }

    // camera(){
    //     position = {0.9f, 0, 0};
    //     look_dir = {-1, 0, 0};
    //     look_dir = linalg::normalize(look_dir);
    //     up = {0, 1, 0};
    //     fov = 90;
    //     aspect = 1;
    // }


    Vec3f geRayDir(float u, float v);
};

class Scene {
public:
    Scene();
    ~Scene();

    static const int LIGHT_MATERIAL_ID = 0;
    static const int BLUE_MATERIAL_ID = 1;
    static const int RED_MATERIAL_ID = 2;
    static const int WHITE_MATERIAL_ID = 3;
    static const int CLOTH_MATERIAL_ID = 4;

    // Add object to scene
    void addObject(obj_data* object, int material_id, bool invert_normals = false);

    // Add an area light
    // p: position
    // n: normal
    // l: side length
    // color: emission vector
    void setAreaLight(Vec3f p, Vec3f n, float l, Vec3f color);
    Vec3f getEmission(Vec3f p, Vec3f inbound_dir);
    bool hitsAreaLight(Vec3f p, Vec3f dir, muni::RayTracer::Octree *octree);
    std::tuple<Vec3f, float, float> sampleAreaLight(Vec3f p);

    ImageMat* renderImage(int w, int h, int spp);
    void renderMaps(int w, int h, ImageMat* &normal_map, ImageMat* &depth_map, ImageMat* &albedo_map);


private: 
    camera cam;

    std::vector<triangle> triangles;

    Vec3f area_light_p;
    Vec3f area_light_n;
    Vec3f area_light_color;
    float area_light_l;
    int area_light_idx;

    std::vector<Material*> materials;

    Vec3f shade_pixel(triangle tri, Vec3f p, Vec3f wo, muni::RayTracer::Octree *octree);
};