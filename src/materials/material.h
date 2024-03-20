#pragma once
#include "../lib/common.h"

const int LIGHT_TYPE = 0;
const int DIFFUSE_TYPE = 1;
const int CLOTH_TYPE = 2;


class Material {
public:
    Material() = default;
    virtual Vec3f eval(const Vec3f& wo, const Vec3f& wi, const Vec3f &n, Vec3f t1 = {0.0f, 1.0f, 0.0f}, Vec3f t2 = {0.0f, 0.0f, 1.0f}) const = 0;
    virtual std::tuple<Vec3f, float> sample(const Vec3f& n) const = 0;
    virtual int type() const = 0;
};