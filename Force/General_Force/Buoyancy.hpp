#pragma once
#include "Force\ExternalForce.hpp"
#include "Grid_Construction/Grid_And_Particle_System.hpp"
/**
 * @class BuoyancyForce
 * @brief 一个实现了Boussinesq近似的浮力模型.
 * 这个力是驱动云层向上生长的核心引擎.
 * 它基于密度和温度与环境的差异来计算垂直方向的加速度.
*/
class BuoyancyForce : public ExternalForce {
private:
    float m_alpha; // 密度系数
    float m_beta;  // 温度系数
    float m_ambientTemperature; // 环境温度
public:
    // 构造函数，初始化浮力参数
    BuoyancyForce(float alpha, float beta, float ambientTemperature);
    // 重写apply函数，应用浮力到MAC网格
    void apply(MACGrid& grid, float dt) override;
};