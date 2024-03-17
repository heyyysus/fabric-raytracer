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
        position = {0, 0, -0.5};
        look_dir = {0, 0, 1};
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

    ImageMat* renderImage(int w, int h);
    ImageMat* renderNormalMap(int w, int h);


private: 
    camera cam;
    std::vector<triangle> triangles;
};