#pragma once

#include "../lib/common.h"
#include "../lib/math_helpers.h"
#include "../lib/sampler.h"
#include "material.h"

class ClothMaterial : public Material {

private:
    Vec3f color;
    float eta = 1.3f;
    // float gamma_s = 12.0f * (M_PI / 180.0f);
    // float gamma_v = 24.0f * (M_PI / 180.0f);
    float gamma_s = 0.1f;
    float gamma_v = 0.5f;
    float kd = 0.1f;

public:

    static Vec3f getTangentVector(float x) {

        return getPitchedVector(x, 10);
    }

    ClothMaterial() : color(Vec3f(0.0f)){ }
    ClothMaterial(Vec3f color) : color(color){}

    Vec3f eval(const Vec3f &wo, const Vec3f&wi, Vec3f t = {0.0f, 1.0f, 0.0f}) const {
        // Implement the light scattering function based on the empirical model

        // Normalize directions
        Vec3f nwo = normalize(wo);
        Vec3f nwi = normalize(wi);
        Vec3f nt = normalize(t);

        // Compute some common vectors and angles
        Vec3f h = normalize(nwi + nwo); // Halfway vector
        float cosThetaI = dot(nwi, nt);
        float cosThetaO = dot(nwo, nt);
        float cosThetaH = dot(h, nt);

        // Calculate Fresnel reflectance (simplified version)
        float Fr = fresnel(cosThetaI, eta);

        // Surface and volume scattering components
        float surfaceScattering = gaussian(cosThetaH, gamma_s) * Fr;
        float volumeScattering = gaussian(cosThetaH, gamma_v) * (1.0f - kd) + kd;

        // Combine surface and volume scattering with the base color
        Vec3f fs = this->color * (surfaceScattering + volumeScattering);

        return fs;
    }

    int type() const {
        return CLOTH_TYPE;
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

private:
    float gaussian(float x, float stddev) const {
        return exp(-(x * x) / (2.0f * stddev * stddev)) / (stddev * sqrt(2.0f * M_PI));
    }

    float fresnel(float cosTheta, float eta) const {
        return pow((eta - 1.0f) / (eta + 1.0f), 2.0f);
    }

    static Vec3f getPitchedVector(float x, int numOscillations) {

        float scaledX = x * 2.0f * M_PI * numOscillations;
    
        float xComponent = std::sin(scaledX);
        float yComponent = std::sqrt(1.0f - xComponent * xComponent);
        float zComponent = 0.0f;

        Vec3f pitchedVector = Vec3f(xComponent, yComponent, zComponent);

        return normalize(pitchedVector);
    }

};