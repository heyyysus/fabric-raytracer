#include "scene.h"
#include "lib/image.h"
#include "lib/ray_tracer.h"
#include <limits>

Scene::Scene(){
    // Initialize camera
    this->cam = camera();
}


ImageMat* Scene::renderNormalMap(int w, int h){
    ImageMat* img = new ImageMat(w, std::vector<Vec3f>(h, Vec3f(0, 0, 0)));
    float aspect_ratio = (float)w / (float)h;
    float scale = tan(this->cam.fov * 0.5 * M_PI / 180);
    float inv_w = 1.0 / w;
    float inv_h = 1.0 / h;

    muni::RayTracer::Octree octree = muni::RayTracer::Octree(this->triangles);

    for (int i = 0; i < w; i++){
        if (i % 10 == 0){
            std::cout << "Rendering row: " << i << "/" << w << std::endl;
        }
        for (int j = 0; j < h; j++){
            float u = (2 * i * inv_w - 1) * aspect_ratio * scale;
            float v = (1 - 2 * j * inv_h) * scale;
            Vec3f dir = this->cam.geRayDir(u, v);
            auto [hit, tri, t] = muni::RayTracer::closest_hit(
                this->cam.position, 
                dir, 
                octree,
                this->triangles
            );
            if (hit){
                Vec3f n = t.n;
                (*img)[i][j] = (n + 1) / 2;
                // (*img)[i][j] = {1.0f, 1.0f, 1.0f};
            }
        }
    }

    return img;
}


Vec3f camera::geRayDir(float u, float v){

    float fov_rad = fov * M_PI / 180;
    float half_height = tan(fov_rad / 2);
    float half_width = aspect * half_height;

    Vec3f right = cross(look_dir, up);
    Vec3f up_ = cross(right, look_dir);

    Vec3f image_center = position + look_dir;
    Vec3f image_up = up_ * half_height;
    Vec3f image_right = right * half_width;

    Vec3f image_point = image_center + image_right * u + image_up * v;
    return linalg::normalize(image_point - position);
}

void Scene::addObject(obj_data* object, int material_id){
    for (auto& face : object->faces){
        triangle tri;
        tri.v0 = object->vertices[std::get<0>(face[0]) - 1];
        tri.v1 = object->vertices[std::get<0>(face[1]) - 1];
        tri.v2 = object->vertices[std::get<0>(face[2]) - 1];
        tri.n = linalg::normalize(linalg::cross(tri.v1 - tri.v0, tri.v2 - tri.v0));
        tri.material_id = material_id;
        this->triangles.push_back(tri);
    }
}