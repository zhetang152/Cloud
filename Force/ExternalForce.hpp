#pragma once
#include "Grid_Construction\Grid_And_Particle_System.hpp"

class ExternalForce {
public:
    virtual ~ExternalForce() = default;
    // apply函数接收grid和时间步长dt作为参数
    virtual void apply(MACGrid& grid, float dt) = 0;
};