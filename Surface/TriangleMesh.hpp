#pragma once
#include "Grid_Construction\GridAndParticleSystem.hpp"
#include <vector>

// 用于存储最终三角网格的数据结构
struct TriangleMesh {
    std::vector<Vector3D> vertices;
    std::vector<Vector3D> normals; //用于平滑着色
    std::vector<int> faces;       //每3个整数代表一个三角形面
};