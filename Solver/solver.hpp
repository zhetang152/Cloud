#pragma once
#include "Grid_Construction/grid.hpp"
#include "Grid_Construction/Grid_And_Particle_System.hpp"
#include <iostream>
#include <cmath>
namespace Solver {
    //基本的向量操作
    //dot product
    float dotProduct(const Grid<float>& a, const Grid<float>& b, const Grid<CellType>&cellTypes);
    //scalar product & vector addition
    void saxpy(Grid<float>& y, float a, const Grid<float>& x, const Grid<CellType>& cellTypes);
    //PCG
    //matrix-vector product
    void applyA(
    Grid<float>& result, 
    const Grid<float>& p,
    const Grid<CellType>& cellTypes,
    const Grid<float>& Adiag,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k
    );
    //计算MACGrid的离散散度
    Grid<float> discrete_divergence(const MACGrid& grid);
    //buildmatrixA
    void buildMatrixA(
        Grid<float>& Adiag, 
        Grid<float>& Aplus_i, 
        Grid<float>& Aplus_j, 
        Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes,
        float dx, 
        float dt, 
        float rho
    );
    float dotProduct_FVM(const Grid<float>& a, const Grid<float>& b, const MACGrid& grid);
    void saxpy_FVM(Grid<float>& y, float a, const Grid<float>& x, const MACGrid& grid);
    //[FVM & Variation] discrete divergence
    Grid<float> discrete_divergence_FVM(const MACGrid& grid);
    //[FVM & Variation] build matrix A
    void buildMatrixA_FVM(
        Grid<float>& Adiag, 
        Grid<float>& Aplus_i, 
        Grid<float>& Aplus_j, 
        Grid<float>& Aplus_k,
        const MACGrid& grid,
        float dt
    );
    void applyA_FVM(
        Grid<float>& result, 
        const Grid<float>& p,
        const MACGrid& grid,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k
    );
    //preconditioner by MIC(0)
    void MIC0preconditioner(
        Grid<float>& precon,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    );
    //apply preconditioner
    void applyPreconditioner(
    Grid<float>& z,
    const Grid<float>& r,
    const Grid<float>& precon,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k,
    const Grid<CellType>& cellTypes
    );
    void MIC0preconditioner_FVM(
    Grid<float>& precon,
    const Grid<float>& Adiag,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k,
    const MACGrid& grid // 传入 MACGrid 以获取 volumeFractions
);
void Solver::applyPreconditioner_FVM(
    Grid<float>& z,
    const Grid<float>& r,
    const Grid<float>& precon,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k,
    const MACGrid& grid
);
    //PCG solver
    struct SystemMatrix {
        Grid<float> Adiag;
        Grid<float> Aplus_i;
        Grid<float> Aplus_j;
        Grid<float> Aplus_k;
        Grid<float> precon;
        SystemMatrix(int nx, int ny, int nz) :
            Adiag(nx, ny, nz, 0.0f),
            Aplus_i(nx, ny, nz, 0.0f),
            Aplus_j(nx, ny, nz, 0.0f),
            Aplus_k(nx, ny, nz, 0.0f),
            precon(nx, ny, nz, 0.0f)
        {}
    };
    void PCG(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        float dx,
        int maxIterations,
        float tolerance
    );
    void PCG_FVM(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const MACGrid& grid,
        int maxIterations,
        float tolerance
    );
}