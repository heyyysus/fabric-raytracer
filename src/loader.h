#pragma once

#include "lib/common.h"
#include "lib/math_helpers.h"

#include <string>
#include <iostream>
#include <fstream>



struct obj_data {
    std::vector<Vec3f> vertices;
    std::vector<Vec2f> tex_coords;
    std::vector<Vec3f> normals;
    std::vector<std::vector<std::tuple<int, int, int> > > faces;
};

struct mtl_data {
    std::string name;
    std::string texture;
};

obj_data* load_obj(const std::string& filePath);

void normalize_obj_vertices(obj_data* data);

void apply_transformation(obj_data* data, const std::vector<std::vector<float> >& transformation);
