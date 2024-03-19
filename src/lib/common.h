#pragma once
#include "linalg.h"

template<int N, class T> using Vec = linalg::vec<T, N>;
template<class T> using Vec3 = Vec<3, T>;
template<class T> using Vec2 = Vec<2, T>;


using Vec3f = Vec3<float>;
using Vec2f = Vec2<float>;

struct triangle {
    Vec3f v0, v1, v2;
    Vec2f uv0, uv1, uv2;
    Vec3f n;
    int material_id;

    triangle(){
        v0 = {0, 0, 0};
        v1 = {0, 0, 0};
        v2 = {0, 0, 0};
        n = {0, 0, 0};
        material_id = -1;
    }

    triangle(std::tuple<Vec3f, Vec3f, Vec3f, Vec3f, int> tri){

        v0 = std::get<0>(tri);
        v1 = std::get<1>(tri);
        v2 = std::get<2>(tri);
        n = std::get<3>(tri);;
        material_id = std::get<4>(tri);
    }

    // static Vec3f get_tangent_vector(triangle &tri, Vec3f &p) {

    //     Vec3f E0 = tri.v1 - tri.v0;
    //     Vec3f E1 = tri.v2 - tri.v0;

    //     Vec2f dUV0 = tri.uv1 - tri.uv0;
    //     Vec2f dUV1 = tri.uv2 - tri.uv0;

    //     float f = 1.0f / (dUV0.x * dUV1.y - dUV1.x * dUV0.y);

    //     Vec3f tangent = f * (dUV1.y * E0 - dUV0.y * E1);

    //     return normalize(tangent);

    // }

    static Vec3f get_tangent_vector(triangle &tri, Vec3f &p) {

        float scaledx = p.z * 4.0f * M_PI;

        float x = sin(scaledx);
        float y = cos(scaledx);
        float z = 0.0f;

        return normalize(Vec3f(x, y, z));

    }

    static std::tuple<bool, float>
    ray_triangle_intersect(triangle tri, Vec3f ray_origin, Vec3f ray_direction,
                           float t_min, float t_max);
    
};
