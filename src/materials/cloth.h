#pragma once

#include "../lib/common.h"
#include "../lib/math_helpers.h"
#include "../lib/sampler.h"
#include "material.h"
#include <cmath>

class ClothMaterial : public Material {

private:
    Vec3f color;
    float eta = 1.46f;
    float gamma_s = 12.0f * (M_PI / 180.0f);
    float gamma_v = 24.0f * (M_PI / 180.0f);
    // float gamma_s = 0.1f;
    // float gamma_v = 0.5f;
    float kd = 0.3f;

public:

    ClothMaterial() : color(Vec3f(0.0f)){ }
    ClothMaterial(Vec3f color) : color(color){}

    Vec3f eval(const Vec3f &wo, const Vec3f&wi, Vec3f t = {0.0f, 1.0f, 0.0f}) const {

        Vec3f nwo = normalize(wo);
        Vec3f nwi = normalize(wi);
        Vec3f nt = normalize(t);

        Vec3f h = normalize(nwi + nwo); 
        float cosThetaI = dot(nwi, nt);
        float cosThetaO = dot(nwo, nt);
        float cosThetaH = dot(h, nt);

        float thetaD = (cosThetaI - cosThetaO)/2;
        float cosThetaD = cosf(thetaD);

        float Fr = fresnel(cosThetaI, eta);

        float surfaceScattering = gaussian(cosThetaH, gamma_s) * Fr;
        float volumeScattering = gaussian(cosThetaH, gamma_v) * (1.0f - kd) + kd;

        Vec3f fs = this->color * INV_PI * (1.0f / (cosThetaD * cosThetaD)) * (surfaceScattering + volumeScattering);

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

        float scaledX = x * 2.0f * M_PI * (float)(numOscillations);
    
        float xComponent = std::sinf(scaledX);
        float yComponent = 0.0f;
        float zComponent = std::cosf(scaledX);

        Vec3f pitchedVector = Vec3f(xComponent, yComponent, zComponent);

        return normalize(pitchedVector);
    }

    std::vector<Vec3f> sampleTangents(const Vec3f& p, const Vec3f& n, const std::vector<double>& offsets, double length) {

        std::vector<Vec3f> tangents;

        Vec3f baseTangent(1, 0, 0);
        
        Vec3f arbitrary(0.0, 1.0, 0.0);
        if (std::fabs(dot(n, arbitrary)) > 0.999) { 
            arbitrary = Vec3f(1.0, 0.0, 0.0);
        }

        Vec3f rotationAxis = cross(n, arbitrary);
        rotationAxis = normalize(rotationAxis);

        for (double offset : offsets) {

            Vec3f rotatedTangent = muni::rotate(baseTangent, rotationAxis, offset);
            rotatedTangent = normalize(rotatedTangent);

            rotatedTangent.x *= length;
            rotatedTangent.y *= length;
            rotatedTangent.z *= length;

            tangents.push_back(rotatedTangent);
        }

        return tangents;
    }

};