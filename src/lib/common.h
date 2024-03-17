#pragma once
#include "linalg.h"

template<int N, class T> using Vec = linalg::vec<T, N>;
template<class T> using Vec3 = Vec<3, T>;
template<class T> using Vec2 = Vec<2, T>;


using Vec3f = Vec3<float>;
using Vec2f = Vec2<float>;

struct triangle {
    Vec3f v0, v1, v2;
    Vec3f n;
    int material_id;

    triangle(){
        v0 = {0, 0, 0};
        v1 = {0, 0, 0};
        v2 = {0, 0, 0};
        n = {0, 0, 0};
        material_id = -1;
    }

    static std::tuple<bool, float>
    ray_triangle_intersect(triangle tri, Vec3f ray_origin, Vec3f ray_direction,
                           float t_min, float t_max);
};
