#include "loader.h"
#include "lib/linalg.h"

obj_data* load_obj(const std::string& filePath) {
    std::ifstream file(filePath);
    obj_data* data = nullptr;

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return data;
    }

    data = new obj_data();

    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == 'v') {
            if (line[1] == ' ') {
                Vec3f vertex;
                sscanf(line.c_str(), "v %f %f %f", &vertex.x, &vertex.y, &vertex.z);
                data->vertices.push_back(vertex);
            } else if (line[1] == 't') {
                Vec2f tex_coord;
                sscanf(line.c_str(), "vt %f %f", &tex_coord.x, &tex_coord.y);
                data->tex_coords.push_back(tex_coord);
            } else if (line[1] == 'n') {
                Vec3f normal;
                sscanf(line.c_str(), "vn %f %f %f", &normal.x, &normal.y, &normal.z);
                data->normals.push_back(normal);
            }
        } else if (line[0] == 'f') {
            std::vector<std::tuple<int, int, int> > f1, f2;
            int v1, v2, v3, v4, t1, t2, t3, t4, n1, n2, n3, n4;
            sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &v1, &t1, &n1, &v2, &t2, &n2, &v3, &t3, &n3, &v4, &t4, &n4);
            
            // Create 2 triangles out of the quad
            f1.push_back(std::make_tuple(v1, t1, n1));
            f1.push_back(std::make_tuple(v2, t2, n2));
            f1.push_back(std::make_tuple(v3, t3, n3));
            data->faces.push_back(f1);

            f2.push_back(std::make_tuple(v1, t1, n1));
            f2.push_back(std::make_tuple(v3, t3, n3));
            f2.push_back(std::make_tuple(v4, t4, n4));
            data->faces.push_back(f2);
        }
    }

    file.close();

    return data;
}

// Average length of vertices to 1
void normalize_obj_vertices(obj_data* data){
    Vec3f min = {INFINITY, INFINITY, INFINITY};
    Vec3f max = {-INFINITY, -INFINITY, -INFINITY};

    for (Vec3f vertex : data->vertices){
        min = linalg::min(min, vertex);
        max = linalg::max(max, vertex);
    }

    Vec3f center = (min + max) / 2;
    float max_len = linalg::length(max - min);

    for (Vec3f& vertex : data->vertices){
        vertex = (vertex - center) / max_len;
    }
}


void apply_transformation(obj_data* data, const std::vector<std::vector<float> >& transformation){
    for (Vec3f& vertex : data->vertices){
        vertex = muni::matrix_multiply(vertex, transformation);
    }
}