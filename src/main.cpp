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
        for (int i = 0; i < width; i++){
            for (int j = 0; j < height; j++){
                for (int k = 0; k < 3; k++){
                    float val = muni::clamp(static_cast<float>(255.0f * this->data[j][i][k]), 0.0f, 255.0f);
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
                Vec3f pixel = tone_map_Aces(data[i % height][i / height]);
                _data[3 * i + j] = static_cast<uint8_t>(
                    // 255.0f * std::max(0.0f, std::min(1.0f, pixels[i][j])));
                    255.0f * std::max(0.0f, std::min(1.0f, pixel[j])));
            }
        }
        // Save the image to a png file.
        return stbi_write_png(filename.c_str(), width, height, 3, _data.data(),
                              sizeof(uint8_t) * 3 * width) != 0;
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
        {{3, 5, 5}, {8, 5, 5}, {4, 5, 5}},
        {{5, 6, 6}, {1, 6, 6}, {4, 6, 6}},
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

obj_data* create_fabric(float scale){
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

        v[0] -= 0.95f;
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

    return fabric;
}

int main(int argc, char** argv){

    int w = 256;
    int h = 256;
    int spp = 16;
    float light_intensity = 70.0f;

    // obj_data* object = load_obj("models/" + fn + ".obj");
    obj_data* object = load_obj("models/jacket.obj");
    obj_data* walls0 = create_walls(1.0f);
    obj_data* walls1 = create_walls(1.0f);

    walls0->faces.erase(walls0->faces.begin(), walls0->faces.end()-2);
    walls1->faces.erase(walls1->faces.end()-2, walls1->faces.end());

    std::vector<std::vector<float> > rot = {
        {0, 1, 0},
        {-1, 0, 0},
        {0, 0, 1},
    };

    apply_transformation(object, rot);

    // Apply 180 degree rotation around the y-axis
    std::vector<std::vector<float> > rot2 = {
        {-1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    apply_transformation(object, rot2);

    normalize_obj_vertices(object);

    // for (Vec3f& v : object->vertices){
    //     v[2] -= 0.5f;
    // }

    if (object == nullptr){
        return 1;
    }

    std::cout << "Loaded " << object->vertices.size() << " vertices" << std::endl;
    std::cout << "Loaded " << object->tex_coords.size() << " texture coordinates" << std::endl;
    std::cout << "Loaded " << object->normals.size() << " normals" << std::endl;


    Scene* scene = new Scene();

    scene->addObject(object, Scene::CLOTH_MATERIAL_ID, true);
    scene->addObject(walls0, Scene::BLUE_MATERIAL_ID);
    scene->addObject(walls1, Scene::WHITE_MATERIAL_ID);

    Vec3f light_pos = {0.95f, 0.0f, 0.0f};
    Vec3f light_dir = { -1.0f, 0.0f, 0.0f };

    scene->setAreaLight(light_pos, light_dir, 0.3f, Vec3f(light_intensity));

    ImageMat *nmap, *dmap, *amap;

    std::cout << "Rendering..." << std::endl;
    scene->renderMaps(w, h, nmap, dmap, amap);

    std::cout << "ImageMat to Image..." << std::endl;
    Image normal_map(*nmap), depth_map(*dmap), albedo_map(*amap);

    std::cout << "Saving..." << std::endl;
    normal_map.save("out/normal.png");
    depth_map.save("out/depth.png");
    albedo_map.save("out/albedo.png");

    Image img(*scene->renderImage(w, h, spp));
    img.save_with_tonemapping("out/origin.png");



    return 0;
}