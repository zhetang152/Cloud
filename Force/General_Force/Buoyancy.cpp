#include "Buoyancy.hpp"
#include <cmath> // For std::pow

BuoyancyForce::BuoyancyForce(float alpha, float beta, float ambientTemperature):
    m_alpha(alpha), m_beta(beta), m_ambientTemperature(ambientTemperature){}

void BuoyancyForce::apply(MACGrid& grid, float dt) {
    int nx = grid.getDimX();
    int ny = grid.getDimY();
    int nz = grid.getDimZ();

    // 遍历所有的垂直速度分量w
    auto& w = grid.w();
    const auto& density = grid.density();
    const auto& temperature = grid.temperature();

    for(int k = 1; k < nz; ++k) {
        for(int j = 0; j < ny; ++j) {
            for(int i = 0; i < nx; ++i) {
                // 计算w(i,j,k)所在位置的平均温度和密度
                float avg_temp = (temperature(i, j, k) + temperature(i, j, k - 1)) * 0.5f;
                float avg_density = (density(i, j, k) + density(i, j, k - 1)) * 0.5f;
                // Boussinesq近似
                float buoyancy = -m_alpha * avg_density + m_beta * (avg_temp - m_ambientTemperature);
                w(i, j, k) += buoyancy * dt;
            }
        }
    }
}