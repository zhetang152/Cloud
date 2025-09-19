#pragma once
#include "Grid_Construction\Grid_And_Particle_System.hpp"
namespace FLIPSolver{
    /**
     *@brief: P2G 
    */
    void ParticleToGrid(MACGrid& grid);
    /**
     *@brief: G2P与FLIP&PIC混合
     *@param grid: MACGrid
     *@param u_old: 上一帧的网格速度
     *@param v_old: 上一帧的网格速度
     *@param w_old: 上一帧的网格速度
     *@param alpha: FLIP与PIC的混合比例
    */
   void GridToParticle(MACGrid& grid, const Grid<float>& u_old, const Grid<float>& v_old, const Grid<float>& w_old, float alpha=0.95f);
}