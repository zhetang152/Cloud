#include "repo.hpp"
float smoothstep(float edge0, float edge1, float x) {

    x = std::max(0.0f, std::min(1.0f, (x - edge0) / (edge1 - edge0)));

    return x * x * (3.0f - 2.0f * x);

}
int main() {
    //1. 模拟参数
    const int resolution = 64; //网格分辨率
    const float dt = 0.01f; //时间步长
    const int totalframes = 200; //总时长
    //物理参数
    const float ambient_temp = 273.15f; //环境温度
    const float buoyancy_alpha = 0.05f; //密度对浮力的影响系数
    const float buoyancy_beta = 0.35f; //温度对浮力的影响系数

    //2. 初始化
    MACGrid grid(resolution, resolution, resolution,1.0f/resolution);
    grid.temperature().fill(ambient_temp);
    //播种一小块湿热空气作为云种子
    int centerX = resolution / 2;
    int centerY = resolution / 2;
    int seed_height = resolution / 8;
    int seed_radius = resolution / 4;
    for(int k = 0; k < resolution; ++k){
        for(int j = 0; j < resolution; ++j){
            for(int i = 0; i < resolution; ++i){
                float dist_sq = pow(i - centerX, 2) + pow(j - centerY, 2);
                if(dist_sq < pow(seed_radius, 2)){
                    grid.density()(i, j, k) = 1.0f;
                    //物理模拟, 将此处温度调高
                    grid.temperature()(i, j, k) = ambient_temp + 5.0f;
                }
            }
        }
    }
    //创建一个包含所有力的列表
    std::vector<std::unique_ptr<ExternalForce>> forces;
    forces.push_back(std::make_unique<BuoyancyForce>(buoyancy_alpha, buoyancy_beta, ambient_temp));
    //3. 主循环
    for(int frame = 0; frame < totalframes; ++frame){
        //a. 应用外力
        std::cout << "施加外力...\n";
        for(const auto& force: forces){
            force-> apply(grid, dt);
        }
        //b. 平流模拟
        std::cout << "进行速度场平流...\n";
        Advector::advect_velocity(grid, dt);
        std::cout << "进行密度场平流...\n";
        grid.density() = Advector::advect(grid.density(), grid, dt);
        std::cout << "进行温度场平流...\n";
        grid.temperature() = Advector::advect(grid.temperature(), grid, dt);

        //c. 投影步骤
        std::cout << "进行投影步骤...\n";
        int nx = grid.getDimX();
        int ny = grid.getDimY();
        int nz = grid.getDimZ();
        float dx = grid.getDx();
        Grid<float> divergence = Solver::discrete_divergence(grid);

        Grid<float>& pressure = grid.pressure();
        pressure.fill(0.0f);
        
        const float rho = 1.0f;
        Solver::SystemMatrix matrix(nx, ny, nz);
        Solver::buildMatrixA(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, 
                             grid.celltypes(), dx, dt, rho);
        Solver::PCG(pressure, divergence, matrix, grid.celltypes(), dx, 200, 1e-5f);

        float scale = dt / rho;
        auto& u = grid.u();
        auto& v = grid.v();
        auto& w = grid.w();
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 1; i < nx; ++i) {
            if (grid.celltypes()(i, j, k) == CellType::FLUID || grid.celltypes()(i - 1, j, k) == CellType::FLUID) {
                u(i, j, k) -= scale * (pressure(i, j, k) - pressure(i - 1, j, k)) / dx;
            }
        }
        for (int k = 0; k < nz; ++k) for (int j = 1; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            if (grid.celltypes()(i, j, k) == CellType::FLUID || grid.celltypes()(i, j - 1, k) == CellType::FLUID) {
                v(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j - 1, k)) / dx;
            }
        }
        for (int k = 1; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
             if (grid.celltypes()(i, j, k) == CellType::FLUID || grid.celltypes()(i, j, k - 1) == CellType::FLUID) {
                w(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j, k - 1)) / dx;
            }
        }

        //d. 输出当前帧数据
        std::cout << "输出当前帧数据...\n";
        MarchingCubes mc(resolution, resolution, resolution);
        mc.init_all();
        for(int k = 0; k < resolution; ++k){
            for(int j = 0; j < resolution; ++j){
                for(int i = 0; i < resolution; ++i){
                    //使用体积分数作为等值面提取依据
                    mc.set_data(grid.density()(i, j, k) - 0.5f, i, j, k);
                }
            }
        }
        mc.run(0.5f);

        if(mc.nverts() > 0){
            char filename[100];
            sprintf(filename, "D:/Code/XLAB/Fluid/result/cloud_frame_%04d.ply", frame);
            mc.writePLY(filename, false);
        }
        mc.clean_all();
    }
    std::cout << "模拟完成！\n";
    return 0;
}