#include"solver.hpp"
#include<vector>

namespace Solver {
    float dotProduct(const Grid<float>& a, const Grid<float>& b, const Grid<CellType>& cellTypes) {
    float result = 0.0f;
    int nx = a.getWidth();
    int ny = a.getHeight();
    int nz = a.getDepth();
    for(int k=0;k<nz;++k){
        for(int j=0;j<ny;++j){
            for(int i=0;i<nx;++i){
                if(cellTypes(i, j, k) == CellType::FLUID) {
                    result += a(i, j, k) * b(i, j, k);
                }
            }
        }
    }
    return result;
    }
    void saxpy(Grid<float>& y, float a, const Grid<float>& x, const Grid<CellType>& cellTypes) {
        int nx = y.getWidth();
        int ny = y.getHeight();
        int nz = y.getDepth();
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i, j, k) == CellType::FLUID) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }
    Grid<float> discrete_divergence(const MACGrid& grid){
        int nx = grid.celltypes().getWidth();
        int ny = grid.celltypes().getHeight();
        int nz = grid.celltypes().getDepth();
        Grid<float> negetivedivergence(nx,ny,nz,0.0f);
        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        const auto& celltypes = grid.celltypes();
        for(int k = 0; k<nz;++k){
            for(int j = 0; j<ny;++j){
                for(int i = 0;i<nx;++i){
                    if(celltypes(i,j,k)==CellType::FLUID){
                        float u_right = u(i+1,j,k);
                        float u_left = u(i,j,k);
                        float v_top = v(i,j+1,k);
                        float v_bottom = v(i,j,k);
                        float w_front = w(i,j,k+1);
                        float w_back = w(i,j,k);
                        float divergence = (u_right-u_left)+(v_top-v_bottom)+(w_front-w_back);
                        negetivedivergence(i,j,k)=-divergence;
                    }
                }
            }
        }
        return negetivedivergence;
    }
    void buildMatrixA(
    Grid<float>& Adiag, 
    Grid<float>& Aplus_i, 
    Grid<float>& Aplus_j, 
    Grid<float>& Aplus_k,
    const Grid<CellType>& cellTypes,
    float dx, 
    float dt, 
    float rho
    ) {
        int nx = Adiag.getWidth();
        int ny = Adiag.getHeight();
        int nz = Adiag.getDepth();
        bool pressure_pinned = false;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    // 初始化，确保每次循环都是干净的
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;

                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        if (!pressure_pinned) {
                            Adiag(i, j, k) = 1.0f;
                            pressure_pinned = true;
                            continue;
                        }

                        float diag_val = 0;
                        if (i < nx - 1 && cellTypes(i + 1, j, k) != CellType::SOLID) {
                            diag_val++;
                            if (cellTypes(i + 1, j, k) == CellType::FLUID) Aplus_i(i, j, k) = -1.0f;
                        }
                        if (i > 0 && cellTypes(i - 1, j, k) != CellType::SOLID) diag_val++;
                        if (j < ny - 1 && cellTypes(i, j + 1, k) != CellType::SOLID) {
                            diag_val++;
                            if (cellTypes(i, j + 1, k) == CellType::FLUID) Aplus_j(i, j, k) = -1.0f;
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) != CellType::SOLID) diag_val++;
                        if (k < nz - 1 && cellTypes(i, j, k + 1) != CellType::SOLID) {
                            diag_val++;
                            if (cellTypes(i, j, k + 1) == CellType::FLUID) Aplus_k(i, j, k) = -1.0f;
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) != CellType::SOLID) diag_val++;
                        
                        Adiag(i, j, k) = diag_val;
                    }
                }
            }
        }
    }
    void applyA(
        Grid<float>& result, 
        const Grid<float>& p,
        const Grid<CellType>& cellTypes,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k
    ) {
        int nx = result.getWidth();
        int ny = result.getHeight();
        int nz = result.getDepth();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        
                        // result_i = Σ (A_ij * p_j)
                        //         = A_ii*p_i + Σ (A_in * p_n) for neighbors n

                        // 1. 对角线部分 A(i,j,k),(i,j,k) * p(i,j,k)
                        float val = Adiag(i, j, k) * p(i, j, k);

                        // 2. 非对角线部分（邻居的影响）
                        // X方向
                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) {
                            val += Aplus_i(i, j, k) * p(i + 1, j, k);
                        }
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            // 利用对称性: A(i,i-1) = A(i-1,i) = Aplus_i(i-1)
                            val += Aplus_i(i - 1, j, k) * p(i - 1, j, k);
                        }

                        // Y方向
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) {
                            val += Aplus_j(i, j, k) * p(i, j + 1, k);
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            val += Aplus_j(i, j - 1, k) * p(i, j - 1, k);
                        }

                        // Z方向
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) {
                            val += Aplus_k(i, j, k) * p(i, j, k + 1);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            val += Aplus_k(i, j, k - 1) * p(i, j, k - 1);
                        }

                        result(i, j, k) = val;

                    } else {
                        result(i, j, k) = 0.0f; // 非流体单元格结果为零
                    }
                }
            }
        }
    }
    float dotProduct_FVM(const Grid<float>& a, const Grid<float>& b, const MACGrid& grid) {
        float result = 0.0f;
        int nx = a.getWidth();
        int ny = a.getHeight();
        int nz = a.getDepth();
        const auto& volumeFractions = grid.volumeFractions();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (volumeFractions(i, j, k) > 0.0f) {
                        result += a(i, j, k) * b(i, j, k) * volumeFractions(i, j, k);
                    }
                }
            }
        }
        return result;
    }
    void saxpy_FVM(Grid<float>& y, float a, const Grid<float>& x, const MACGrid& grid){
        int nx = y.getWidth();
        int ny = y.getHeight();
        int nz = y.getDepth();
        const auto& volumeFractions = grid.volumeFractions();
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    // ✅ 判断标准改为体积分数
                    if(volumeFractions(i, j, k) > 0.0f) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }
    Grid<float> discrete_divergence_FVM(const MACGrid& grid){
        int nx = grid.getDimX();
        int ny = grid.getDimY();
        int nz = grid.getDimZ();
        float dx = grid.getDx();
        Grid<float> negetivedivergence(nx, ny, nz, 0.0f);
        
        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        
        const auto& volumeFractions = grid.volumeFractions();
        const auto& solidVelocity = grid.solidvelocity();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (volumeFractions(i, j, k) > 0.0f) {
                        float fluid_flux = 
                            (u(i + 1, j, k) * u_area(i + 1, j, k) - u(i, j, k) * u_area(i, j, k)) +
                            (v(i, j + 1, k) * v_area(i, j + 1, k) - v(i, j, k) * v_area(i, j, k)) +
                            (w(i, j, k + 1) * w_area(i, j, k + 1) - w(i, j, k) * w_area(i, j, k));
                        float solid_flux = 0.0f;
                        // X方向
                        //进行边界检查
                        int i_p1 = std::min(i + 1, nx - 1);
                        int i_m1 = std::max(i - 1, 0);
                        Vector3D solid_vel_right = 0.5f * (solidVelocity(i_p1, j, k) + solidVelocity(i, j, k));
                        Vector3D solid_vel_left  = 0.5f * (solidVelocity(i,j,k) + solidVelocity(i_m1,j,k));
                        solid_flux += solid_vel_right.x * (1.0f - u_area(i + 1, j, k));
                        solid_flux -= solid_vel_left.x  * (1.0f - u_area(i, j, k));
                        // Y方向
                        int j_p1 = std::min(j + 1, ny - 1);
                        int j_m1 = std::max(j - 1, 0);
                        Vector3D solid_vel_top    = 0.5f * (solidVelocity(i, j_p1, k) + solidVelocity(i, j, k));
                        Vector3D solid_vel_bottom = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j_m1, k));
                        solid_flux += solid_vel_top.y    * (1.0f - v_area(i, j + 1, k));
                        solid_flux -= solid_vel_bottom.y * (1.0f - v_area(i, j, k));

                        // Z方向
                        int k_p1 = std::min(k + 1, nz - 1);
                        int k_m1 = std::max(k - 1, 0);
                        Vector3D solid_vel_front = 0.5f * (solidVelocity(i, j, k_p1) + solidVelocity(i, j, k));
                        Vector3D solid_vel_back = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j, k_m1));
                        solid_flux += solid_vel_front.z * (1.0f - w_area(i, j, k + 1));
                        solid_flux -= solid_vel_back.z * (1.0f - w_area(i, j, k));
                        // 计算负散度
                        negetivedivergence(i, j, k) = 
                            -(fluid_flux - solid_flux) / dx;
                    }
                }
            }
        }
        return negetivedivergence;
    };
    void buildMatrixA_FVM(
        Grid<float>& Adiag, 
        Grid<float>& Aplus_i, 
        Grid<float>& Aplus_j, 
        Grid<float>& Aplus_k,
        const MACGrid& grid,
        float dt
    ){
        int nx = grid.getDimX();
        int ny = grid.getDimY();
        int nz = grid.getDimZ();
        float dx = grid.getDx();
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        const auto& volume_frac = grid.volumeFractions();
        const auto& density = grid.density();
        for (int k = 0; k < nz; ++k){
            for (int j = 0; j < ny; ++j){
                for (int i = 0; i < nx; ++i){
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;
                    if (volume_frac(i ,j, k) > 0.0f){
                        float diag_val = 0.0f;
                        //x方向
                        // 右侧
                        if (i < nx -1){
                            float rho_face = 0.5f * (density(i, j, k) + density(i + 1, j, k));
                            float scale_face = dt / (rho_face * dx * dx);
                            Aplus_i(i, j, k) = -scale_face * u_area(i + 1, j, k);
                            diag_val += scale_face * u_area(i + 1, j, k);
                        }
                        // 左侧
                        if (i > 0){
                            float rho_face = 0.5f * (density(i, j, k) + density(i - 1, j, k));
                            float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * u_area(i, j, k);
                        }
                        //y方向
                        // 上侧
                        if (j < ny - 1) {
                            float rho_face = 0.5f * (density(i, j, k) + density(i, j + 1, k));
                            float scale_face = dt / (rho_face * dx * dx);
                            Aplus_j(i, j, k) = -scale_face * v_area(i, j + 1, k);
                            diag_val += scale_face * v_area(i, j + 1, k);
                        }
                        // 下侧
                        if (j > 0) {
                            float rho_face = 0.5f * (density(i, j, k) + density(i, j - 1, k));
                            float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * v_area(i, j, k);
                        }
                        //z方向
                        // 前侧
                        if (k < nz - 1) {
                            float rho_face = 0.5f * (density(i, j, k) + density(i, j, k + 1));
                            float scale_face = dt / (rho_face * dx * dx);
                            Aplus_k(i, j, k) = -scale_face * w_area(i, j, k + 1);
                            diag_val += scale_face * w_area(i, j, k + 1);
                        }
                        // 后侧
                        if (k > 0) {
                            float rho_face = 0.5f * (density(i, j, k) + density(i, j, k - 1));
                            float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * w_area(i, j, k);
                        }
                        Adiag(i, j, k) = diag_val;
                    }
                }
            }
        };
    }
    void applyA_FVM(
        Grid<float>& result, 
        const Grid<float>& p,
        const MACGrid& grid,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k
    ) {
        int nx = result.getWidth();
        int ny = result.getHeight();
        int nz = result.getDepth();
        const auto& volumeFraction = grid.volumeFractions();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (grid.volumeFractions()(i, j, k) > 0.0f) {
                        // 1. 对角线部分 A(i,j,k),(i,j,k) * p(i,j,k)
                        float val = Adiag(i, j, k) * p(i, j, k);

                        // 2. 非对角线部分（邻居的影响）
                        // X方向
                        if (i < nx - 1 && volumeFraction(i + 1, j, k) > 0.0f) {
                        val += Aplus_i(i, j, k) * p(i + 1, j, k);
                        }
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            val += Aplus_i(i - 1, j, k) * p(i - 1, j, k);
                        }

                        // Y方向
                        if (j < ny - 1 && volumeFraction(i, j + 1, k) > 0.0f) {
                            val += Aplus_j(i, j, k) * p(i, j + 1, k);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            val += Aplus_j(i, j - 1, k) * p(i, j - 1, k);
                        }

                        // Z方向
                        if (k < nz - 1 && volumeFraction(i, j, k + 1) > 0.0f) {
                            val += Aplus_k(i, j, k) * p(i, j, k + 1);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            val += Aplus_k(i, j, k - 1) * p(i, j, k - 1);
                        }

                        result(i, j, k) = val;

                    } else {
                        result(i, j, k) = 0.0f; // 非流体单元格结果为零
                    }
                }
            }
        }
    }
    void MIC0preconditioner(
        Grid<float>& precon,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    ){
        int nx = precon.getWidth();
        int ny = precon.getHeight();
        int nz = precon.getDepth();
        const float tau = 0.97f; // MIC(0)的调节参数
        const float sgm = 0.25f; // 安全参数
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        //IC(0)
                        float term_i_sq = 0.0f, term_j_sq = 0.0f, term_k_sq = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            float precon_val = precon(i - 1, j, k);
                            term_i_sq = (Aplus_i(i - 1, j, k) * precon_val) * (Aplus_i(i - 1, j, k) * precon_val);
                        }
                        if (j> 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            float precon_val = precon(i, j - 1, k);
                            term_j_sq = (Aplus_j(i, j - 1, k) * precon_val) * (Aplus_j(i, j - 1, k) * precon_val);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            float precon_val = precon(i, j, k - 1);
                            term_k_sq = (Aplus_k(i, j, k - 1) * precon_val) * (Aplus_k(i, j, k - 1) * precon_val);
                        }
                        //MIC(0)
                        float comp_i = 0.0f, comp_j = 0.0f, comp_k = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            float precon_sq = precon(i - 1, j, k) * precon(i - 1, j, k);
                            comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            float precon_sq = precon(i, j - 1, k) * precon(i, j - 1, k);
                            comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            float precon_sq = precon(i, j, k - 1) * precon(i, j, k - 1);
                            comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                        }
                        float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
                        if (e < sgm * Adiag(i, j, k)) {
                        e = Adiag(i, j, k);
                        }
                        // 存储e的平方根的倒数
                        precon(i, j, k) = 1.0f / std::sqrt(e);
                    } else {
                        precon(i, j, k) = 0.0f;
                    }
                }
            }
        }
    }
    void applyPreconditioner(
        Grid<float>& z,
        const Grid<float>& r,
        const Grid<float>& precon,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    ){
        int nx = z.getWidth();
        int ny = z.getHeight();
        int nz = z.getDepth();
        Grid<float> q(nx, ny, nz, 0.0f);
        //前向替换 (Forward Substitution), 求解 Lq = r
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        float offdiag_sum = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                        }
                        float t = r(i, j, k) - offdiag_sum;
                        q(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
        //后向替换 (Backward Substitution), 求解 L^Tz = q
        for(int k=nz-1;k>=0;--k){
            for(int j=ny-1;j>=0;--j){
                for(int i=nx-1;i>=0;--i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        float offdiag_sum = 0.0f;
                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                        }
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                        }
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                        }
                        float t = q(i, j, k) - offdiag_sum;
                        z(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
    }
    void MIC0preconditioner_FVM(
        Grid<float>& precon,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const MACGrid& grid // 传入 MACGrid 以获取 volumeFractions
    ){
        int nx = precon.getWidth();
        int ny = precon.getHeight();
        int nz = precon.getDepth();
        const float tau = 0.97f;
        const float sgm = 0.25f;

        const auto& volumeFraction = grid.volumeFractions();

        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    //判断标准改为体积分数
                    if(volumeFraction(i, j, k) > 0.0f){
                        
                        float term_i_sq = 0.0f, term_j_sq = 0.0f, term_k_sq = 0.0f;

                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            float precon_val = precon(i - 1, j, k);
                            term_i_sq = (Aplus_i(i - 1, j, k) * precon_val) * (Aplus_i(i - 1, j, k) * precon_val);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            float precon_val = precon(i, j - 1, k);
                            term_j_sq = (Aplus_j(i, j - 1, k) * precon_val) * (Aplus_j(i, j - 1, k) * precon_val);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            float precon_val = precon(i, j, k - 1);
                            term_k_sq = (Aplus_k(i, j, k - 1) * precon_val) * (Aplus_k(i, j, k - 1) * precon_val);
                        }
                        
                        float comp_i = 0.0f, comp_j = 0.0f, comp_k = 0.0f;
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            float precon_sq = precon(i - 1, j, k) * precon(i - 1, j, k);
                            comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            float precon_sq = precon(i, j - 1, k) * precon(i, j - 1, k);
                            comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            float precon_sq = precon(i, j, k - 1) * precon(i, j, k - 1);
                            comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                        }

                        float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
                        
                        if (e < sgm * Adiag(i, j, k)) {
                            e = Adiag(i, j, k);
                        }
                        
                        precon(i, j, k) = 1.0f / std::sqrt(e);

                    } else {
                        precon(i, j, k) = 0.0f;
                    }
                }
            }
        }
    }
    void applyPreconditioner_FVM(
        Grid<float>& z,
        const Grid<float>& r,
        const Grid<float>& precon,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const MACGrid& grid
    ){
        int nx = z.getWidth();
        int ny = z.getHeight();
        int nz = z.getDepth();
        Grid<float> q(nx, ny, nz, 0.0f);
        
        const auto& volumeFraction = grid.volumeFractions();

        // 前向替换, 求解 Lq = r
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    //判断标准改为体积分数
                    if(volumeFraction(i, j, k) > 0.0f){
                        float offdiag_sum = 0.0f;
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                        }
                        float t = r(i, j, k) - offdiag_sum;
                        q(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }

        // 后向替换, 求解 L^Tz = q
        for(int k = nz - 1; k >= 0; --k){
            for(int j = ny - 1; j >= 0; --j){
                for(int i = nx - 1; i >= 0; --i){
                    if(volumeFraction(i, j, k) > 0.0f){
                        float offdiag_sum = 0.0f;
                        if (i < nx - 1 && volumeFraction(i + 1, j, k) > 0.0f) {
                            offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                        }
                        if (j < ny - 1 && volumeFraction(i, j + 1, k) > 0.0f) {
                            offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                        }
                        if (k < nz - 1 && volumeFraction(i, j, k + 1) > 0.0f) {
                            offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                        }
                        float t = q(i, j, k) - offdiag_sum;
                        z(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
    }
    void PCG(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        float dx,
        int maxIterations,
        float tolerance
    ){
        //获取网格尺寸
        int nx = p.getWidth();
        int ny = p.getHeight();
        int nz = p.getDepth();
        //初始化
        Grid<float> r(nx, ny, nz, 0.0f);
        Grid<float> Ap(nx, ny, nz, 0.0f);
        r = b;
        //计算初始残差
        float initial_residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
        if (initial_residual_norm < 1e-9) {
            std::cout << "Initial residual is zero."<<std::endl;
            return;
        }
        //z=M-1 * r
        Grid<float> z(nx, ny, nz, 0.0f);
        applyPreconditioner(z,r,matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
        //d = z
        Grid<float> d = z;
        // delta_new = r · z
        float delta_new = dotProduct(r, z, cellTypes);

        //PCG
        for(int k =0; k < maxIterations;++k){
            //q = A * d
            Grid<float> q(nx, ny, nz, 0.0f);
            applyA(q, d, cellTypes, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
            //alpha = delta_new / (d · q)
            float d_q = dotProduct(d, q, cellTypes);
            float alpha = delta_new / d_q;
            //p = p + alpha * d
            saxpy(p, alpha, d, cellTypes);
            //r = r - alpha * q
            saxpy(r, -alpha, q, cellTypes);
            //收敛性检查
            float residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
            std::cout << "Iteration " << k + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
                return;
            }
            //z_new = M-1 * r_new
            applyPreconditioner(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
            float delta_old = delta_new;
            delta_new = dotProduct(r, z, cellTypes);
            //beta = delta_new / delta_old
            float beta = delta_new / delta_old;
            //d = z + beta * d
            Grid<float> temp_d = d;
            d = z;
            saxpy(d, beta, temp_d, cellTypes);
        }
        std::cout << "PCG did not converge after" << maxIterations << std::endl;
    };
    void PCG_FVM(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const MACGrid& grid,
        int maxIterations,
        float tolerance
    ){
        int nx = p.getWidth();
        int ny = p.getHeight();
        int nz = p.getDepth();
        // 初始化残差
        Grid<float> r(nx, ny, nz, 0.0f);
        r = b;
        // 计算初始残差范数用于收敛判断
        //调用 FVM 版本的dotProduct
        float initial_residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
        if (initial_residual_norm < 1e-9) {
            std::cout << "FVM: Initial residual is zero." << std::endl;
            return;
        }
        // z = M^-1 * r
        Grid<float> z(nx, ny, nz, 0.0f);
        //调用 FVM 版本的applyPreconditioner
        applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
        // d = z
        Grid<float> d = z;
        // delta_new = r · z
        float delta_new = dotProduct_FVM(r, z, grid);
        // --- 主循环 ---
        for (int k = 0; k < maxIterations; ++k) {
            // q = A * d
            Grid<float> q(nx, ny, nz, 0.0f);
            applyA_FVM(q, d, grid, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
            
            // alpha = delta_new / (d · q)
            float d_q = dotProduct_FVM(d, q, grid);
            float alpha = delta_new / d_q;

            // p = p + alpha * d
            saxpy_FVM(p, alpha, d, grid);

            // r = r - alpha * q
            saxpy_FVM(r, -alpha, q, grid);

            // 收敛性检查
            float residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
            std::cout << "FVM Iteration " << k + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "FVM Converged after " << k + 1 << " iterations." << std::endl;
                return;
            }
            
            // z_new = M^-1 * r_new
            applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
            
            float delta_old = delta_new;
            //调用 FVM 版本的 dotProduct
            delta_new = dotProduct_FVM(r, z, grid);

            // beta = delta_new / delta_old
            float beta = delta_new / delta_old;
            
            // d = z + beta * d
            Grid<float> temp_d = d;
            d = z;
            saxpy_FVM(d, beta, temp_d, grid);
        }
        
        std::cout << "FVM PCG did not converge after " << maxIterations << " iterations." << std::endl;
    }
}