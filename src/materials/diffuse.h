#pragma once

#include "../lib/common.h"

class DiffuseMaterial {

private:
    Vec3f color;

public:

    DiffuseMaterial() : color(Vec3f(0.0f)){}
    DiffuseMaterial(Vec3f color) : color(color){}

    Vec3f eval() const {
        return this->color;
    }

};