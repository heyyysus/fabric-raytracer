#include "scene.h"
#include "lib/common.h"
#include "lib/image.h"
#include "lib/linalg.h"
#include "lib/math_helpers.h"
#include "lib/ray_tracer.h"
#include "materials/diffuse.h"
#include "materials/cloth.h"
#include <limits>
#include <tuple>
#include "lib/sampler.h"

#define EPS 1e-6f

Scene::Scene(camera cam){
    // Initialize camera
    this->cam = cam;

    this->area_light_idx = -1;

    Material *blue_material, *red_material, *white_material, *cloth_material;
    blue_material = new DiffuseMaterial(Vec3f(0.1f, 0.2f, 0.5f));
    red_material = new DiffuseMaterial(Vec3f(0.85f, 0, 0));
    white_material = new DiffuseMaterial(Vec3f(0.7f, 0.7f, 0.7f));
    cloth_material = new ClothMaterial(Vec3f(0.3f, 0.05f, 0.1f), Vec3f(0.1f, 0.05f, 0.3f));

    this->materials.push_back(blue_material);
    this->materials.push_back(red_material);
    this->materials.push_back(white_material);
    this->materials.push_back(cloth_material);

    muni::UniformSampler::init(0);
}

std::tuple<bool, float> triangle::ray_triangle_intersect(triangle tri, Vec3f ray_origin, Vec3f ray_direction,
        float t_min, float t_max){

    const Vec3f abs_ray_direction = abs(ray_direction);
        unsigned int axis = 0;
        if (abs_ray_direction[1] > abs_ray_direction[0] &&
            abs_ray_direction[1] > abs_ray_direction[2])
            axis = 1;
        if (abs_ray_direction[2] > abs_ray_direction[0] &&
            abs_ray_direction[2] > abs_ray_direction[1])
            axis = 2;

        unsigned int kz = axis;
        unsigned int kx = (kz + 1) % 3;
        unsigned int ky = (kx + 1) % 3;
        if (ray_direction[kz] < 0.0f) {
            unsigned int swap = kx;
            kx = ky;
            ky = swap;
        }

        float Sx = ray_direction[kx] / ray_direction[kz];
        float Sy = ray_direction[ky] / ray_direction[kz];
        float Sz = 1.f / ray_direction[kz];

        const Vec3f A = tri.v0 - ray_origin;
        const Vec3f B = tri.v1 - ray_origin;
        const Vec3f C = tri.v2 - ray_origin;

        const float Ax = A[kx] - Sx * A[kz];
        const float Ay = A[ky] - Sy * A[kz];
        const float Bx = B[kx] - Sx * B[kz];
        const float By = B[ky] - Sy * B[kz];
        const float Cx = C[kx] - Sx * C[kz];
        const float Cy = C[ky] - Sy * C[kz];

        float U = Cx * By - Cy * Bx;
        float V = Ax * Cy - Ay * Cx;
        float W = Bx * Ay - By * Ax;

        if (U == 0.f || V == 0.f || W == 0.f) {
            double CxBy = static_cast<double>(Cx) * static_cast<double>(By);
            double CyBx = static_cast<double>(Cy) * static_cast<double>(Bx);
            U = (float)(CxBy - CyBx);
            double AxCy = static_cast<double>(Ax) * static_cast<double>(Cy);
            double AyCx = static_cast<double>(Ay) * static_cast<double>(Cx);
            V = (float)(AxCy - AyCx);
            double BxAy = static_cast<double>(Bx) * static_cast<double>(Ay);
            double ByAx = static_cast<double>(By) * static_cast<double>(Ax);
            W = (float)(BxAy - ByAx);
        }

        if ((U < 0.f || V < 0.f || W < 0.f) && (U > 0.f || V > 0.f || W > 0.f))
            return {false, 0.0f};

        float det = U + V + W;
        if (det == 0.f) return {false, 0.0f};

        const float Az = Sz * A[kz];
        const float Bz = Sz * B[kz];
        const float Cz = Sz * C[kz];
        const float T = U * Az + V * Bz + W * Cz;
        const float rcp_det = 1.f / det;

        const float t = T * rcp_det;
        const Vec3f barycentrics = Vec3f(U, V, W) * rcp_det;

        if (t < t_min || t > t_max) return {false, 0.0f};

        return {true, t};

}


void Scene::renderMaps(int w, int h, ImageMat* &normal, ImageMat* &depth, ImageMat* &albedo){
    normal = new ImageMat(w, std::vector<Vec3f>(h, Vec3f(0, 0, 0)));
    depth = new ImageMat(w, std::vector<Vec3f>(h, Vec3f(0, 0, 0)));
    albedo = new ImageMat(w, std::vector<Vec3f>(h, Vec3f(0, 0, 0)));

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
            auto [hit, t_min, tri] = muni::RayTracer::closest_hit(
                this->cam.position, 
                dir, 
                octree,
                this->triangles
            );
            if (hit){
                if (tri.material_id == 0){
                    (*normal)[i][j] = {1.0f, 1.0f, 1.0f};
                    (*depth)[i][j] = linalg::clamp(Vec3f(1.0f - t_min), Vec3f(0.0f), Vec3f(1.0f));
                    (*albedo)[i][j] = this->area_light_color;
                } else {    
                    Vec3f n = tri.n;

                    Vec3f p = this->cam.position + t_min * dir;

                    Vec3f t_v = Vec3f(0.0f);
                    if (tri.material_id == Scene::CLOTH_MATERIAL_ID){
                        ClothMaterial* mat = (ClothMaterial*)this->materials.at(tri.material_id - 1);
                        t_v = mat->sample_tangent(tri.n);
                    }

                    (*normal)[i][j] = (n + 1) / 2;
                    (*depth)[i][j] = linalg::clamp(Vec3f(1.0f - t_min), Vec3f(0.0f), Vec3f(1.0f));
                    (*albedo)[i][j] = this->materials.at(tri.material_id - 1)->eval(-dir, normalize(Vec3f{0.0f, -0.5f, -0.5f}), tri.n, t_v) * M_PI;
                    // if (tri.material_id == Scene::CLOTH_MATERIAL_ID)
                    //     (*albedo)[i][j] = linalg::clamp(t_v, Vec3f(0.0f), Vec3f(1.0f));
                }
            }
        }
    }
}

Vec3f Scene::shade_pixel(triangle tri, Vec3f p, Vec3f wo, muni::RayTracer::Octree *octree){
    Vec3f pixel = {0, 0, 0};

    Vec3f L_dir = Vec3f(0.0f);

    auto [w_light, pdf_light, l_dist] = this->sampleAreaLight(p);

    float t_max = l_dist - 2.0f * EPS;

    bool lightBlocked = 
        muni::RayTracer::any_hit(p + EPS * tri.n, w_light, t_max, *octree, this->triangles);


    Vec3f t_v = Vec3f(0.0f);
    Vec3f normal = tri.n;
    Vec3f eval_d = Vec3f(0.0f);
    Vec3f eval_i = Vec3f(0.0f);

    Vec3f wi;
    float pdf;

    if(tri.material_id != Scene::LIGHT_MATERIAL_ID){
        auto [wi_s, pdf_s] = this->materials.at(tri.material_id - 1)->sample(tri.n);
        wi = wi_s;
        pdf = pdf_s;
    }

    if (tri.material_id == Scene::CLOTH_MATERIAL_ID){
        bool switch_t = false;
        ClothMaterial* mat = (ClothMaterial*)this->materials.at(tri.material_id - 1);
        if (normal != Vec3f(1.0f, 0.0f, 0.0f)){
            normal = linalg::normalize(Vec3f({ 0.0f, p.y, p.z }));
        }
        t_v = mat->sample_tangent(normal);
        Vec3f t_v_2 = normalize(cross(t_v, normal));

        // if(t_v_2 

        // if (abs(length(t_v_2) - 0.5f) > EPS) 
        //     std::cout << "Tangent vector is not orthogonal to normal" << std::endl;

        eval_d = mat->eval(-wo, w_light, normal, t_v, t_v_2);
        eval_i = mat->eval(-wo, wi, normal, t_v, t_v_2);
    } else if (tri.material_id != Scene::LIGHT_MATERIAL_ID){
        
        eval_d = this->materials.at(tri.material_id - 1)->eval(wo, w_light, normal, t_v);
        eval_i = this->materials.at(tri.material_id - 1)->eval(wo, wi, normal, t_v);

    }

    float cos_theta_light_normal = linalg::dot(normal, w_light);

    if (!lightBlocked && cos_theta_light_normal >= 0){

        Vec3f eval = this->getEmission(p, w_light) * cos_theta_light_normal / pdf_light;

        eval *= (tri.material_id == 0) ? Vec3f(1.0f) : eval_d;

        L_dir += eval;

    }

    Vec3f L_indir = Vec3f(0.0f);

    float ksi = muni::UniformSampler::next1d();
    // float ksi = 1.0f;

    if(ksi <= 0.8f && tri.material_id != 0){

        Material* mat = this->materials.at(tri.material_id - 1);
        auto [hit, t_min, hit_triangle] = muni::RayTracer::closest_hit(p + EPS * tri.n, wi, *octree, this->triangles);
        Vec3f q = p + t_min * wi;

        if(hit && mat->type() != LIGHT_TYPE){

            Vec3f shade = eval_i;

            shade *= shade_pixel(hit_triangle, q, -wi, octree) * linalg::dot(normal, wi);
            shade *= (1.0f/pdf) * (1.0f / 0.8f);
            L_indir += shade;
        }
    }

    return L_dir + L_indir;
}


ImageMat* Scene::renderImage(int w, int h, int spp){
    ImageMat* img = new ImageMat(w, std::vector<Vec3f>(h, Vec3f(0, 0, 0)));
    float aspect_ratio = (float)w / (float)h;
    float scale = tan(this->cam.fov * 0.5 * M_PI / 180);
    float inv_w = 1.0 / w;
    float inv_h = 1.0 / h;


    muni::RayTracer::Octree octree = muni::RayTracer::Octree(this->triangles);

    for (int i = 0; i < w; i++){
        if (i % 10 == 0){
            std::cout << "Rendering col: " << i << "/" << w << std::endl;
        }
        for (int j = 0; j < h; j++){
            float u = (2 * (i + muni::UniformSampler::next1d()) * inv_w - 1) * aspect_ratio * scale;
            float v = (1 - 2 * (j + muni::UniformSampler::next1d()) * inv_h) * scale;

            

            Vec3f dir = this->cam.geRayDir(u, v);

            for (int k = 0; k < spp; k++){
                auto [hit, t_min, tri] = muni::RayTracer::closest_hit(
                    this->cam.position, 
                    dir, 
                    octree,
                    this->triangles
                );
                if (hit){
                    Vec3f p = this->cam.position + t_min * dir;
                    Vec3f n = tri.n;

                    if (tri.material_id == 0){
                        (*img)[i][j] += this->getEmission(p, dir);
                    }
                    else {
                        Vec3f pixel = this->shade_pixel(tri, p, dir, &octree);
                        (*img)[i][j] += linalg::clamp(pixel, Vec3f(0.0f), Vec3f(100.0f));
                    }
                } 
            }

            (*img)[i][j] /= (float)(spp);

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

void Scene::addObject(obj_data* object, int material_id, bool invert_normals){
    for (auto& face : object->faces){
        triangle tri;

        tri.v0 = object->vertices[std::get<0>(face[0]) - 1];
        tri.v1 = object->vertices[std::get<0>(face[1]) - 1];
        tri.v2 = object->vertices[std::get<0>(face[2]) - 1];

        tri.n = linalg::normalize(linalg::cross(tri.v1 - tri.v0, tri.v2 - tri.v0));

        if (invert_normals){
            tri.n = -tri.n;
        }

        tri.material_id = material_id;

        tri.uv0 = object->tex_coords[std::get<1>(face[0]) - 1];
        tri.uv1 = object->tex_coords[std::get<1>(face[1]) - 1];
        tri.uv2 = object->tex_coords[std::get<1>(face[2]) - 1];

        this->triangles.push_back(tri);
    }
}

Vec3f Scene::getEmission(Vec3f p, Vec3f inbound_dir){
    if (linalg::dot(inbound_dir, this->area_light_n) < 0){
        return this->area_light_color;
    } else {
        return {0, 0, 0};
    }
}

bool Scene::hitsAreaLight(Vec3f p, Vec3f dir, muni::RayTracer::Octree *octree){

    auto [hit, t_min, tri] = muni::RayTracer::closest_hit(
        p, 
        dir, 
        *octree,
        this->triangles
    );

    return hit && tri.material_id == 0;
}

void Scene::setAreaLight(Vec3f p, Vec3f n, float l, Vec3f color){
    this->area_light_p = p;
    this->area_light_n = n;
    this->area_light_l = l;
    this->area_light_color = color;

    if (this->area_light_idx == -1){
        this->area_light_idx = this->triangles.size();
        this->triangles.push_back(triangle());
        this->triangles.push_back(triangle());
    }

    Vec3f arbitraryVec = (fabs(n.y) < 0.999) ? Vec3f(0, 1, 0) : Vec3f(1, 0, 0); 
    Vec3f u = normalize(cross(n, arbitraryVec)); 
    Vec3f v = normalize(cross(n, u)); 

    Vec3f p0 = p - l * u - l * v;
    Vec3f p1 = p + l * u - l * v;
    Vec3f p2 = p + l * u + l * v;
    Vec3f p3 = p - l * u + l * v;


    triangle t1, t2;
    t1 = std::make_tuple(p0, p1, p2, n, 0);
    t2 = std::make_tuple(p0, p2, p3, n, 0);

    this->triangles.at(this->area_light_idx) = t1;
    this->triangles.at(this->area_light_idx + 1) = t2;
}

std::tuple<Vec3f, float, float> Scene::sampleAreaLight(Vec3f p){
    auto [eps0, eps1] = muni::UniformSampler::next2d();
    Vec3f arbitraryVec = (fabs(area_light_n.y) < 0.999) ? Vec3f(0, 1, 0) : Vec3f(1, 0, 0); 
    Vec3f u = normalize(cross(area_light_n, arbitraryVec)); 
    Vec3f v = normalize(cross(area_light_n, u)); 

    Vec3f sample_point = this->area_light_p + (eps0 - 0.5f) * this->area_light_l * u + (eps1 - 0.5f) * this->area_light_l * v;
    float pdf_area = 1.0f / (this->area_light_l * this->area_light_l);

    Vec3f dir = sample_point - p;
    float dist_sqr = linalg::length2(dir);
    dir = linalg::normalize(dir);

    float pdf_solid_angle = dist_sqr * pdf_area * (1.0f / linalg::dot(-dir, this->area_light_n));

    return {dir, pdf_solid_angle, sqrtf(dist_sqr)};
}

