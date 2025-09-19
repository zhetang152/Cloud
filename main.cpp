#include "repo.hpp"
#include <vector>
#include <memory>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <random> // 用于生成随机云团位置

// fBm (分形布朗运动) 函数，用于生成更自然的噪声
float fBm(const PerlinNoise& noise, float x, float y, float z, int octaves, float persistence) {
    float total = 0.0f;
    float frequency = 1.0f;
    float amplitude = 1.0f;
    float maxValue = 0.0f;

    for(int i = 0; i < octaves; i++) {
        total += noise.getValue(x * frequency, y * frequency, z * frequency) * amplitude;
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= 2.0f;
    }
    return total / maxValue;
}

// smoothstep 函数，用于平滑过渡
float smoothstep(float edge0, float edge1, float x) {
    x = std::max(0.0f, std::min(1.0f, (x - edge0) / (edge1 - edge0)));
    return x * x * (3.0f - 2.0f * x);
}

int main() {
    // 1. 模拟参数
    const int resolution = 128;
    const float dt = 0.01f;
    const int totalframes = 300;
    
    // 物理参数
    const float ambient_temp = 273.15f;
    const float buoyancy_alpha = 0.05f;
    const float buoyancy_beta = 0.35f;

    // 2. 初始化
    MACGrid grid(resolution, resolution, resolution, 1.0f / resolution);

    // 设置边界条件 - Z=0平面为固体地面，其余为开放天空
    for(int k = 0; k < resolution; ++k) {
        for(int j = 0; j < resolution; ++j) {
            for(int i = 0; i < resolution; ++i) {
                if (k == 0) { // Z=0 是地面
                    grid.celltypes()(i, j, k) = CellType::SOLID;
                } else {
                    grid.celltypes()(i, j, k) = CellType::AIR;
                }
            }
        }
    }
    grid.temperature().fill(ambient_temp);
    
    // 在远离地面的高空生成初始云场
    PerlinNoise noise(12345); 

    // --- 噪声参数 ---
    int shape_octaves = 4;
    float shape_persistence = 0.5f;
    int detail_octaves = 6;
    float detail_persistence = 0.5f;

    // --- 云层宏观参数 ---
    int base_height = resolution / 6;      // 云底的平坦高度 (Lifting Condensation Level)
    int max_cloud_height = resolution / 2;
    float temp_amplitude = 35.0f;          // 稍稍增强浮力
    
    // --- “雕刻”参数 ---
    float density_threshold = 0.45f;
    float erosion_strength = 0.75f;

    for (int k = base_height; k < max_cloud_height; ++k) {
        for (int j = 0; j < resolution; ++j) {
            for (int i = 0; i < resolution; ++i) {
                
                // 1. 用低频噪声构建云的大致形状
                float base_shape_noise = fBm(noise, 
                                             (float)i/resolution * 2.5f, 
                                             (float)j/resolution * 2.5f, 
                                             (float)k/resolution * 3.5f, 
                                             shape_octaves, shape_persistence);
                base_shape_noise = (base_shape_noise + 1.0f) * 0.5f;

                // 2. MODIFICATION: 创建一个从下到上逐渐减弱的梯度
                float height_gradient = 1.0f - smoothstep((float)base_height, (float)max_cloud_height, (float)k);
                float final_shape = base_shape_noise * height_gradient;

                if (final_shape > density_threshold) {
                    
                    // 3. 用高频噪声来“侵蚀”云的形状，创造细节
                    float detail_noise = fBm(noise, 
                                             (float)i/resolution * 12.0f, 
                                             (float)j/resolution * 12.0f, 
                                             (float)k/resolution * 16.0f, 
                                             detail_octaves, detail_persistence);
                    detail_noise = (detail_noise + 1.0f) * 0.5f;

                    float eroded_density = final_shape - detail_noise * erosion_strength;

                    if (eroded_density > 0) {
                        float intensity = smoothstep(0.0f, 1.0f, eroded_density);
                        grid.celltypes()(i, j, k) = CellType::FLUID;
                        grid.density()(i, j, k) = intensity;
                        grid.temperature()(i, j, k) = ambient_temp + intensity * temp_amplitude * height_gradient; // 温度也受梯度影响
                    }
                }
            }
        }
    }
    
    // 创建力列表
    std::vector<std::unique_ptr<ExternalForce>> forces;
    forces.push_back(std::make_unique<BuoyancyForce>(buoyancy_alpha, buoyancy_beta, ambient_temp));

    // 3. 主循环
    for(int frame = 0; frame < totalframes; ++frame){
        std::cout << "\n--- Frame " << frame << " ---\n";

        // a. 应用外力
        std::cout << "Applying external forces...\n";
        for(const auto& force: forces){
            force->apply(grid, dt);
        }

        // b. 平流
        std::cout << "Advecting fields...\n";
        Advector::advect_velocity(grid, dt);
        grid.density() = Advector::advect(grid.density(), grid, dt);
        grid.temperature() = Advector::advect(grid.temperature(), grid, dt);
        
        // 尺寸变量
        int nx = grid.getDimX(); int ny = grid.getDimY(); int nz = grid.getDimZ();
        float dx = grid.getDx();

        // 速度限制 (CFL)
        float max_vel_sq = 0.0f;
        auto& u = grid.u(); auto& v = grid.v(); auto& w = grid.w();
        for(int k=0; k<nz; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx+1; ++i) max_vel_sq = std::max(max_vel_sq, u(i,j,k)*u(i,j,k));
        for(int k=0; k<nz; ++k) for(int j=0; j<ny+1; ++j) for(int i=0; i<nx; ++i) max_vel_sq = std::max(max_vel_sq, v(i,j,k)*v(i,j,k));
        for(int k=0; k<nz+1; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx; ++i) max_vel_sq = std::max(max_vel_sq, w(i,j,k)*w(i,j,k));
        float max_vel = std::sqrt(max_vel_sq);
        float cfl_limit = dx * 5.0f / dt; 
        if (max_vel > cfl_limit) {
            std::cout << "WARNING: CFL violated (max_vel=" << max_vel << "). Clamping velocity.\n";
            float clamp_scale = cfl_limit / max_vel;
            for(int k=0; k<nz; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx+1; ++i) u(i,j,k) *= clamp_scale;
            for(int k=0; k<nz; ++k) for(int j=0; j<ny+1; ++j) for(int i=0; i<nx; ++i) v(i,j,k) *= clamp_scale;
            for(int k=0; k<nz+1; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx; ++i) w(i,j,k) *= clamp_scale;
        }

        // c. 投影
        std::cout << "Projection step...\n";
        Grid<float> divergence = Solver::discrete_divergence(grid);
        float div_scale = dx / dt;
        for(int k=0; k<nz; ++k) for(int j=0; j<ny; ++j) for(int i=0; i<nx; ++i) {
            divergence(i,j,k) *= div_scale;
        }
        Grid<float>& pressure = grid.pressure();
        pressure.fill(0.0f);
        // 钉死第一个流体单元的散度为0，以确保PCG收敛
        bool pinned = false;
        for (int k = 0; k < nz && !pinned; ++k) {
            for (int j = 0; j < ny && !pinned; ++j) {
                for (int i = 0; i < nx && !pinned; ++i) {
                    if (grid.celltypes()(i, j, k) == CellType::FLUID) {
                        divergence(i, j, k) = 0.0f;
                        pinned = true;
                    }
                }
            }
        }

        const float rho = 1.0f;
        Solver::SystemMatrix matrix(nx, ny, nz);
        
        Solver::buildMatrixA(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes(), dx, dt, rho);
        
        Solver::MIC0preconditioner(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes());
        Solver::PCG(pressure, divergence, matrix, grid.celltypes(), dx, 200, 1e-5f);
        //应用物理尺度来更新速度
        float pressure_scale = dt / (rho * dx); 
        for (int k=0; k<nz; ++k) for (int j=0; j<ny; ++j) for (int i=1; i<nx; ++i) {
            if (grid.celltypes()(i-1, j, k) != CellType::SOLID || grid.celltypes()(i, j, k) != CellType::SOLID)
                u(i,j,k) -= pressure_scale * (pressure(i,j,k) - pressure(i-1,j,k));
        }
        for (int k=0; k<nz; ++k) for (int j=1; j<ny; ++j) for (int i=0; i<nx; ++i) {
            if (grid.celltypes()(i, j-1, k) != CellType::SOLID || grid.celltypes()(i, j, k) != CellType::SOLID)
                v(i,j,k) -= pressure_scale * (pressure(i,j,k) - pressure(i,j-1,k));
        }
        for (int k=1; k<nz; ++k) for (int j=0; j<ny; ++j) for (int i=0; i<nx; ++i) {
             if (grid.celltypes()(i, j, k-1) != CellType::SOLID || grid.celltypes()(i, j, k) != CellType::SOLID)
                w(i,j,k) -= pressure_scale * (pressure(i,j,k) - pressure(i,j,k-1));
        }

        // d. 输出
        std::cout << "Exporting frame data...\n";
        MarchingCubes mc(resolution, resolution, resolution);
        mc.init_all();
        for(int k=0; k<resolution; ++k) for(int j=0; j<resolution; ++j) for(int i=0; i<resolution; ++i) {
            mc.set_data(grid.density()(i, j, k), i, j, k);
        }
        mc.run(0.1f);

        if(mc.nverts() > 0){
            char filename[256];
            sprintf(filename, "D:/Code/XLAB/Fluid/result/cloud_frame_%04d.ply", frame);
            mc.writePLY(filename, false);
        }
        mc.clean_all();
    }
    std::cout << "Simulation complete!\n";
    return 0;
}