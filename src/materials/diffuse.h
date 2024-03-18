#pragma once

#include "../lib/common.h"
#include "../lib/math_helpers.h"
#include "../lib/sampler.h"
#include "material.h"

class DiffuseMaterial : public Material {

private:
    Vec3f color;

public:

    DiffuseMaterial() : color(Vec3f(0.0f)){ }
    DiffuseMaterial(Vec3f color) : color(color){}

    Vec3f eval(const Vec3f &wo, const Vec3f&wi, Vec3f t = {0.0f, 1.0f, 0.0f}) const {
        return this->color / M_PI;
    }

    int type() const {
        return DIFFUSE_TYPE;
    }

    std::tuple<Vec3f, float> sample(const Vec3f& n) const {
        
        // =============================================================================================

        Vec2f u = muni::UniformSampler::next2d();

        float theta = u.x * 2.0f * M_PI;
        float phi = acos(sqrt(u.y)); // Cos

        float x = sin(phi) * cos(theta);
        float y = sin(phi) * sin(theta);
        float z = cos(phi); 

        Vec3f wi = Vec3f{x, y, z};
        wi = muni::from_local(wi, n);
        wi = normalize(wi);

        float pdf = cos(phi) / M_PI; //Cos

        return std::make_tuple(wi, pdf);

        // =============================================================================================
    }
};