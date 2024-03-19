#include "lib/common.h"
#include "lib/math_helpers.h"
#include "loader.h"
#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

#include <iostream>
#include <vector>
#include <string>

#include "scene.h"

struct Image {

    int width;
    int height;
    int channels;
    ImageMat data;

    Image(ImageMat& data){
        this->data = data;
        this->width = data.size();
        this->height = data[0].size();
        this->channels = 3;
    }

    void save(const std::string& filePath) const {
        int width = this->data.size();
        int height = this->data[0].size();
        std::vector<unsigned char> raw_data;
        raw_data.reserve(width * height * 3);
        for (int j = 0; j < height; j++){
            for (int i = 0; i < width; i++){
                for (int k = 0; k < 3; k++){
                    float val = muni::clamp(static_cast<float>(255.0f * this->data[i][j][k]), 0.0f, 255.0f);
                    raw_data.push_back((unsigned char)(val));
                }
            }
        }

        stbi_write_png(filePath.c_str(), width, height, this->channels, raw_data.data(), 0);
    }

    Vec3f tone_map_Aces(const Vec3f value) const {
        Vec3f color = 0.6f * value;
        float A = 2.51;
        float B = 0.03;
        float C = 2.43;
        float D = 0.59;
        float E = 0.14;
        color = (color * (A * color + B)) / (color * (C * color + D) + E);
        return color;
    }

    bool save_with_tonemapping(const std::string &filename) const {
        // Convert the floating-point data to 8-bit per channel.
        std::vector<uint8_t> _data(3 * width * height);
        for (int i = 0; i < width * height; i++) {
            for (int j = 0; j < 3; j++) {
                Vec3f pixel = tone_map_Aces(data[i % width][i / width]);
                _data[3 * i + j] = static_cast<uint8_t>(
                    // 255.0f * std::max(0.0f, std::min(1.0f, pixels[i][j])));
                    255.0f * std::max(0.0f, std::min(1.0f, pixel[j])));
            }
        }
        // Save the image to a png file.
        return stbi_write_png(filename.c_str(), width, height, 3, _data.data(),
                              sizeof(uint8_t) * 3 * width) != 0;
    }

    bool save_hdr(const std::string &filename) const {
        // Convert the floating-point data to 8-bit per channel.
        std::vector<float> _data(3 * width * height);
        for (int i = 0; i < width * height; i++) {
            for (int j = 0; j < 3; j++) {
                Vec3f pixel = data[i % width][i / width];
                _data[3 * i + j] = static_cast<float>(pixel[j]);
                // _data[3 * i + j] = static_cast<float>(
                //     std::max(0.0f, std::min(1.0f, pixel[j])));
            }
        }
        // Save the image to a hdr file.
        return stbi_write_hdr(filename.c_str(), width, height, 3, _data.data());
    }

};

obj_data* create_walls(float scale){
    obj_data* walls = new obj_data();
    walls->vertices = {
        {-1, -1, -1},
        {1, -1, -1},
        {1, 1, -1},
        {-1, 1, -1},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, 1},
        {-1, 1, 1}
    };

    for (auto& v : walls->vertices){
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;
    }

    walls->normals = {
        {0, 0, -1},
        {0, 0, 1},
        {1, 0, 0},
        {-1, 0, 0},
        {0, 1, 0},
        {0, -1, 0}
    };

    walls->faces = {
        {{1, 1, 1}, {2, 1, 1}, {3, 1, 1}},
        {{1, 1, 1}, {3, 1, 1}, {4, 1, 1}},
        {{5, 2, 2}, {8, 2, 2}, {7, 2, 2}},
        {{5, 2, 2}, {7, 2, 2}, {6, 2, 2}},
        {{1, 3, 3}, {5, 3, 3}, {6, 3, 3}},
        {{1, 3, 3}, {6, 3, 3}, {2, 3, 3}},
        {{2, 4, 4}, {6, 4, 4}, {7, 4, 4}},
        {{2, 4, 4}, {7, 4, 4}, {3, 4, 4}},
        {{3, 5, 5}, {7, 5, 5}, {8, 5, 5}},
        {{5, 6, 6}, {1, 6, 6}, {4, 6, 6}},
        {{3, 5, 5}, {8, 5, 5}, {4, 5, 5}},
        {{5, 6, 6}, {4, 6, 6}, {8, 6, 6}}
    };

    walls->tex_coords = {
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1},
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}
    };
    return walls;
}

obj_data* create_fabric(float scale, Vec3f offset){
    obj_data* fabric = new obj_data();

    fabric->vertices = {
        {0, -1, -0.75},
        {0, -1, 0.75},
        {0, 1, 0.75},
        {0, 1, -0.75}
    };

    for (auto& v : fabric->vertices){
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;

        // v[0] -= 0.95f;
    }

    fabric->normals = {
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1}
    };

    fabric->tex_coords = {
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}
    };

    fabric->faces = {
        {{1, 1, 1}, {3, 1, 1}, {2, 1, 1}},
        {{1, 1, 1}, {4, 1, 1}, {3, 1, 1}}
    };

    for (auto& v : fabric->vertices){
        v[0] += offset[0];
        v[1] += offset[1];
        v[2] += offset[2];
    }

    return fabric;
}

obj_data* create_cylinder(float scale, Vec3f offset){
    obj_data* cylinder = new obj_data();

    int n = 32;
    float r = 0.5f;
    float h = 1.0f;

    for (int i = 0; i < n; i++){
        float theta = 2.0f * M_PI * (float)(i) / (float)(n);
        float theta2 = 2.0f * M_PI * (float)(i + 1) / (float)(n);

        Vec3f v0 = {0.0f, r * cos(theta), r * sin(theta)};
        Vec3f v1 = {0.0f, r * cos(theta2), r * sin(theta2)};
        Vec3f v2 = {h / 2.0f, r * cos(theta2), r * sin(theta2)};
        Vec3f v3 = {h / 2.0f, r * cos(theta), r * sin(theta)};

        v0[0] *= scale;
        v0[1] *= scale;
        v0[2] *= scale;

        v1[0] *= scale;
        v1[1] *= scale;
        v1[2] *= scale;

        v2[0] *= scale;
        v2[1] *= scale;
        v2[2] *= scale;

        v3[0] *= scale;
        v3[1] *= scale;
        v3[2] *= scale;

        v0[0] += offset[0];
        v0[1] += offset[1];
        v0[2] += offset[2];

        v1[0] += offset[0];
        v1[1] += offset[1];
        v1[2] += offset[2];

        v2[0] += offset[0];
        v2[1] += offset[1];
        v2[2] += offset[2];

        v3[0] += offset[0];
        v3[1] += offset[1];
        v3[2] += offset[2];

        Vec3f n0 = {cos(theta), 0, sin(theta)};
        Vec3f n1 = {cos(theta2), 0, sin(theta2)};
        Vec3f n2 = {cos(theta2), 0, sin(theta2)};
        Vec3f n3 = {cos(theta), 0, sin(theta)};

        n0 = normalize(n0);
        n1 = normalize(n1);
        n2 = normalize(n2);
        n3 = normalize(n3);

        cylinder->vertices.push_back(v0);
        cylinder->vertices.push_back(v1);
        cylinder->vertices.push_back(v2);
        cylinder->vertices.push_back(v3);
        
        cylinder->normals.push_back(n0);
        cylinder->normals.push_back(n1);
        cylinder->normals.push_back(n2);
        cylinder->normals.push_back(n3);

        cylinder->tex_coords.push_back({(float)(i) / (float)(n), 0});
        cylinder->tex_coords.push_back({(float)(i + 1) / (float)(n), 0});
        cylinder->tex_coords.push_back({(float)(i + 1) / (float)(n), 1});
        cylinder->tex_coords.push_back({(float)(i) / (float)(n), 1});

        int idx = i * 4 + 1;
        cylinder->faces.push_back({{idx + 1, idx + 1, idx + 1}, {idx + 2, idx + 3, idx + 3}, {idx + 3, idx + 2, idx + 2}});
        cylinder->faces.push_back({{idx + 1, idx + 1, idx + 1}, {idx + 3, idx + 2, idx + 2}, {idx + 0, idx + 0, idx + 0}});

    }

    Vec3f topCenter = Vec3f(h / 2.0f, 0.0f, 0.0f) * scale + offset;
    cylinder->vertices.push_back(topCenter);
    int topCenterIdx = cylinder->vertices.size(); // Index of the top center vertex
    cylinder->normals.push_back(Vec3f(1, 0, 0)); // Assuming top faces +x direction
    cylinder->tex_coords.push_back(Vec2f(0.5, 0.5)); // Center of UV for the cap

    for (int i = 0; i < n; i++) {
        float theta = 2.0f * M_PI * i / n;
        Vec3f vert = Vec3f(h / 2.0f, r * cos(theta), r * sin(theta)) * scale + offset;
        
        cylinder->vertices.push_back(vert);
        cylinder->normals.push_back(Vec3f(1, 0, 0)); // Top faces +x direction
        float u = (cos(theta) + 1.0f) * 0.5f; // Map to [0, 1] for UV
        float v = (sin(theta) + 1.0f) * 0.5f; // Map to [0, 1] for UV
        cylinder->tex_coords.push_back(Vec2f(u, v));
    }

    for (int i = 0; i < n; i++) {
        int nextIdx = (i + 1) % n;
        int v1 = topCenterIdx + i + 1; // +1 because topCenterIdx is the first
        int v2 = topCenterIdx + nextIdx + 1;

        cylinder->faces.push_back({{topCenterIdx, topCenterIdx, topCenterIdx}, {v1, v1, v1}, {v2, v2, v2}});
    }


    return cylinder;
}

int main(int argc, char** argv){

    int w = 256;
    int h = 256;
    int spp = 8;
    float light_intensity = 80.0f;

    camera cam;
    cam.position = {-0.5, -0.75, -0.3};
    cam.look_dir = normalize(Vec3f(-1.0f,0.0f, 0.0f) - cam.position);
    cam.up = {1, 0, 0};
    cam.fov = 70;
    // cam.aspect = (float)(w) / (float)(h);
    cam.aspect = 1;

    // obj_data* object = load_obj("models/jacket.obj");
    // obj_data* object = create_fabric(0.5f, {-0.95f, 0.0f, 0.0f});
    obj_data* object = create_cylinder(0.5f, {-1.0f, 0.0f, 0.0f});
    obj_data* walls0 = create_walls(1.0f);
    obj_data* walls1 = create_walls(1.0f);

    walls0->faces.erase(walls0->faces.begin(), walls0->faces.end()-2);
    walls1->faces.erase(walls1->faces.end()-2, walls1->faces.end());

    // std::vector<std::vector<float> > rot = {
    //     {0, 1, 0},
    //     {-1, 0, 0},
    //     {0, 0, 1},
    // };

    // apply_transformation(object, rot);

    // // Apply 180 degree rotation around the y-axis
    // std::vector<std::vector<float> > rot2 = {
    //     {-1, 0, 0},
    //     {0, 1, 0},
    //     {0, 0, 1},
    // };

    // apply_transformation(object, rot2);

    // normalize_obj_vertices(object);

    // for (Vec3f& v : object->vertices){
    //     v[1] += 0.9f;
    // }

    if (object == nullptr){
        return 1;
    }

    std::cout << "Loaded " << object->vertices.size() << " vertices" << std::endl;
    std::cout << "Loaded " << object->tex_coords.size() << " texture coordinates" << std::endl;
    std::cout << "Loaded " << object->normals.size() << " normals" << std::endl;


    Scene* scene = new Scene(cam);

    scene->addObject(object, Scene::CLOTH_MATERIAL_ID);
    scene->addObject(walls0, Scene::BLUE_MATERIAL_ID);
    scene->addObject(walls1, Scene::WHITE_MATERIAL_ID);

    Vec3f light_pos = {0.0f, -0.5f, -0.5f};
    Vec3f light_dir = Vec3f{ -1.0f, 0.0f, 0.0f } - light_pos;

    scene->setAreaLight(light_pos, light_dir, 0.3f, Vec3f(light_intensity));

    ImageMat *nmap, *dmap, *amap;

    std::cout << "Rendering..." << std::endl;
    scene->renderMaps(w, h, nmap, dmap, amap);

    std::cout << "ImageMat to Image..." << std::endl;
    Image normal_map(*nmap), depth_map(*dmap), albedo_map(*amap);

    std::cout << "Saving ldr..." << std::endl;
    normal_map.save("out/normal.png");
    depth_map.save("out/depth.png");
    albedo_map.save("out/albedo.png");

    std::cout << "Saving hdr..." << std::endl;
    normal_map.save("out/normal.hdr");
    depth_map.save("out/depth.hdr");
    albedo_map.save("out/albedo.hdr");

    Image img(*scene->renderImage(w, h, spp));
    img.save_with_tonemapping("out/origin.png");
    img.save("out/origin.hdr");


    return 0;
}