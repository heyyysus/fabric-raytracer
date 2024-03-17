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
                    raw_data.push_back((unsigned char)(255 * this->data[i][j][k]));
                }
            }
        }

        stbi_write_png(filePath.c_str(), width, height, this->channels, raw_data.data(), 0);
    }

};

obj_data* create_walls(){
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

int main(int argc, char** argv){

    if (argc != 2){
        std::cerr << "Usage: " << argv[0] << " <filename> (w/o extension, output to <filename>.png)" << std::endl;
        return 1;
    }

    std::string fn = argv[1];

    obj_data* object = load_obj("models/" + fn + ".obj");
    obj_data* walls = create_walls();


    // Apply 90 degree rotation around the z-axis
    // std::vector<std::vector<float> > rot = {
    //     {0, -1, 0, 0},
    //     {1, 0, 0, 0},
    //     {0, 0, 1, 0},
    //     {0, 0, 0, 1}
    // };

    // Inv of rot 
    std::vector<std::vector<float> > rot = {
        {0, 1, 0, 0},
        {-1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };

    apply_transformation(object, rot);

    // Apply 180 degree rotation around the y-axis
    std::vector<std::vector<float> > rot2 = {
        {-1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, -1, 0},
        {0, 0, 0, 1}
    };

    apply_transformation(object, rot2);

    normalize_obj_vertices(object);

    if (object == nullptr){
        return 1;
    }

    std::cout << "Loaded " << object->vertices.size() << " vertices" << std::endl;
    std::cout << "Loaded " << object->tex_coords.size() << " texture coordinates" << std::endl;
    std::cout << "Loaded " << object->normals.size() << " normals" << std::endl;

    const int material_id = 0;

    Scene* scene = new Scene();
    scene->addObject(object, material_id);
    scene->addObject(walls, material_id);

    Image img(*scene->renderNormalMap(512, 512));
    img.save(fn + "_normal_map.png");


    return 0;
}