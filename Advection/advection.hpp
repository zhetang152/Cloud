#pragma once
#include "Grid_Construction\Grid_And_Particle_System.hpp"
#include "Grid_Construction\grid.hpp"
namespace Advector {
    /**
     * @brief 半Lagrange平流方程
     * @param q_old 旧的标量场
     * @param velocitygrid 速度网格
     * @param dt 时间步长
     * @return 新的标量场 
     */
    Vector3D get_velocity_at(const MACGrid& grid, const Vector3D& pos);
    Grid<float> advect(const Grid<float>& q_old, const MACGrid& velocityGrid, float dt);
    void advect_velocity(MACGrid& grid, float dt);
}