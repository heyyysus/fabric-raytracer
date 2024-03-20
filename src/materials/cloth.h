#pragma once

#include "../lib/common.h"
#include "../lib/math_helpers.h"
#include "../lib/sampler.h"
#include "material.h"
#include <cmath>
#include <iostream>

class ClothMaterial : public Material {

private:
    Vec3f color1;
    Vec3f color2;
    float eta = 1.46f;
    float gamma_s = 12.0f * (M_PI / 180.0f);
    float gamma_v = 24.0f * (M_PI / 180.0f);
    // float gamma_s = 0.1f;
    // float gamma_v = 0.5f;
    float kd = 0.3f;
    std::vector<Vec3f> tangents;

public:


    ClothMaterial() : color1(Vec3f(0.0f)), color2(Vec3f(0.0f)){

        const std::vector<float> tangent_lengths = {1.0f, 1.0f};
        const std::vector<double> tangent_offsets = {-25.0, 25.0};
        compute_tangents(Vec3f(0, 0, 1), tangent_lengths, tangent_offsets);

     }
    ClothMaterial(Vec3f color1, Vec3f color2) : color1(color1), color2(color2){

        const std::vector<float> tangent_lengths = {1.0f, 1.0f};
        const std::vector<double> tangent_offsets = {-25.0, 25.0};
        compute_tangents(Vec3f(0, 0, 1), tangent_lengths, tangent_offsets);

    }

    Vec3f eval(const Vec3f &wo, const Vec3f &wi, const Vec3f &n, Vec3f t1 = {0.0f, 1.0f, 0.0f}, Vec3f t2 = {0.0f, 0.0f, 1.0f}) const {

        float a1 = 0.7f;
        float a2 = 0.3f;

        Vec3f fs1 = f_s(wo, wi, t1, color1, n);
        Vec3f fs2 = f_s(wo, wi, t2, color2, n);

        return fs1 * a1 + fs2 * a2;
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

    Vec3f sample_tangent(const Vec3f &n){
        float u = muni::UniformSampler::next1d();
        int idx = (int)(u * (float)(tangents.size() - 1));
        Vec3f tangent = tangents[idx];
        tangent = muni::from_local(tangent, n);
        // if (dot(tangent, n) > EPS || dot(tangent, n) < -EPS)
        //     std::cout << "Tangent dot: " << dot(tangent, n) << std::endl;
        return tangent;
    }

private:

    Vec3f f_s(const Vec3f &wo, const Vec3f&wi, Vec3f t, Vec3f color, Vec3f n) const {
        Vec3f nwo = normalize(wo);
        Vec3f nwi = normalize(wi);
        Vec3f nt = normalize(t);;

        Vec3f h = normalize(nwi + nwo); 
        float cosThetaI = abs(dot(nwi, nt));
        float cosThetaO = abs(dot(nwo, nt));
        float cosThetaH = (cosThetaI + cosThetaO) / 2.0f;
        float cosThetaD = (cosThetaI - cosThetaO) / 2.0f;

        float phi_r = get_phi(nt, nwo, n);
        float phi_i = get_phi(nt, nwi, n);
        float phi_d = phi_i - phi_r;

        // float Fr_wi = Fresnel_r(eta, nwi, n);
        float Fr_wi = fresnel(cosThetaI, eta);
        float Ft_wi = 1.0f - Fr_wi;
        float Ft_wo = 1.0f - fresnel(cosThetaO, eta);
        float Ft = Ft_wi * Ft_wo;

        float surfaceScattering = gaussian(cosThetaH, gamma_s) * Fr_wi * cosf(phi_d / 2.0f);
        Vec3f volumeScattering = color;
        volumeScattering = (volumeScattering * gaussian(cosThetaH, gamma_v) * (1.0f - kd)) + kd;
        volumeScattering *= 1.0f / (cosThetaI + cosThetaO);
        volumeScattering *= Ft;
        // float volumeScattering = gaussian(cosThetaH, gamma_v) * (1.0f - kd) + kd;

        Vec3f fs = color * INV_PI * (surfaceScattering + volumeScattering);

        for (float &f : fs) {
            if (f < 0.0f) {
                f = 0.0f;
            }
            if (std::isnan(f) || std::isinf(f) || f > 1.0f) {
                f = 1.0f;
            }
        }

        return fs;
    }

    float gaussian(float x, float stddev) const {
        return exp(-(x * x) / (2.0f * stddev * stddev)) / (stddev * sqrt(2.0f * M_PI));
    }

    float fresnel(float cosTheta, float eta) const {
        return pow((eta - 1.0f) / (eta + 1.0f), 2.0f);
    }

    float Fresnel_r(float eta, const Vec3f& w, const Vec3f& n) const {

        Vec3f w_norm = normalize(w);
        Vec3f n_norm = normalize(n);

        float cos_theta = std::max(0.0f, std::min(1.0f, dot(w_norm, n_norm)));
        float R0 = std::pow((eta - 1.0f) / (eta + 1.0f), 2.0f);
        float reflectance = R0 + (1.0f - R0) * std::pow(1.0f - cos_theta, 5.0f);

        return reflectance;
}

    float M(const Vec3f &t, const Vec3f &w_r, const Vec3f &phi_r) const{

    }

    float get_phi(const Vec3f &t, const Vec3f &w_r, const Vec3f &n) const {
        Vec3f nw_r = normalize(w_r);
        Vec3f nt = normalize(t);

        Vec3f w_r_proj = nw_r - dot(nw_r, n) * n;
        w_r_proj = normalize(w_r_proj);

        Vec3f b = cross(n, nt);
        b = normalize(b);

        float cos_phi_r = dot(w_r_proj, t);
        float sin_phi_r = dot(w_r_proj, b); 
        float phi_r = atan2(sin_phi_r, cos_phi_r);

        return phi_r;
    }

    static Vec3f getPitchedVector(float x, int numOscillations) {

        float scaledX = x * 2.0f * M_PI * (float)(numOscillations);
    
        float xComponent = std::sinf(scaledX);
        float yComponent = 0.0f;
        float zComponent = std::cosf(scaledX);

        Vec3f pitchedVector = Vec3f(xComponent, yComponent, zComponent);

        return normalize(pitchedVector);
    }

    void compute_tangents(const Vec3f& n, const std::vector<float>& tangent_lengths, const std::vector<double>& tangent_offsets) {

        std::vector<Vec3f> tangents;

        Vec3f baseTangent(1, 0, 0);
        
        Vec3f arbitrary(0.0, 1.0, 0.0);
        if (std::fabs(dot(n, arbitrary)) > 0.999) { 
            arbitrary = Vec3f(1.0, 0.0, 0.0);
        }

        Vec3f rotationAxis = cross(n, arbitrary);
        rotationAxis = normalize(rotationAxis);

        for (double offset : tangent_offsets) {
            for (float length : tangent_lengths) {

                Vec3f rotatedTangent = muni::rotate(baseTangent, rotationAxis, offset);
                rotatedTangent = normalize(rotatedTangent);

                rotatedTangent.x *= length;
                rotatedTangent.y *= length;
                rotatedTangent.z *= length;

                tangents.push_back(rotatedTangent);

            }
        }

        this->tangents = tangents;

    }

};