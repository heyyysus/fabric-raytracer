#pragma once

#include "lib/image.h"
#include "loader.h"

struct camera {
    Vec3f position;
    Vec3f look_dir;
    Vec3f up;
    float fov;
    float aspect;

    camera(){
        position = {0.25, -0.4, -0.5};
        look_dir = {-0.5, 0.65, 1};
        look_dir = linalg::normalize(look_dir);
        up = {0, 1, 0};
        fov = 90;
        aspect = 1;
    }

    Vec3f geRayDir(float u, float v);
};

class Scene {
public:
    Scene();
    ~Scene();

    // Add object to scene
    void addObject(obj_data* object, int material_id);

    // Add an area light
    // p: position
    // n: normal
    // l: side length
    // color: emission vector
    void setAreaLight(Vec3f p, Vec3f n, float l, Vec3f color);
    Vec3f getEmission(Vec3f p, Vec3f inbound_dir);
    bool hitsAreaLight(Vec3f p, Vec3f dir);

    ImageMat* renderImage(int w, int h);
    ImageMat* renderNormalMap(int w, int h);


private: 
    camera cam;
    std::vector<triangle> triangles;
    Vec3f area_light_p;
    Vec3f area_light_n;
    float area_light_l;
    Vec3f area_light_color;
};