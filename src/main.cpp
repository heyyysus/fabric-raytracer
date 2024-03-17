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

int main(int argc, char** argv){

    if (argc != 4){
        std::cerr << "Usage: " << argv[0] << " <input.obj> <input.mtl> <output.png>" << std::endl;
        return 1;
    }

    std::string obj_path = argv[1];
    std::string mtl_path = argv[2];
    std::string output_path = argv[3];

    obj_data* object = load_obj(obj_path);

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

    Image img(*scene->renderNormalMap(512, 512));
    img.save(output_path);


    return 0;
}