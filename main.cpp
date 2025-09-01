#include<iostream>
#include<chrono>
#include<vector>
#include <fstream> // 用于文件读写
#include "NOISE\PerlinNoise.hpp"
#include "Surface\MarchingCubes.h"

int main() {
    // 1. 使用种子创建Perlin噪声生成器
    unsigned int seed = 1338; // 使用一个固定种子，方便每次生成同样的结果
    // unsigned int seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    PerlinNoise perlin(seed);
    std::cout << "使用种子: " << seed << "\n";

    // 2. 定义生成的参数
    const int resolution = 64; // 在一个 64x64x64 的立方体中采样
    const float frequency = 5.0f;    // 噪声频率，控制云朵特征的大小
    const float threshold = 0.5f;   // 密度阈值, 作为marching cubes的iso值

    //3. 初始化MarchingCubes对象
    MarchingCubes mc(resolution, resolution, resolution);
    mc.init_all();

    std::cout << "正在生成OBJ文件，分辨率 " << resolution << "^3 ...\n";

    // 4. 遍历采样区域的每一个点
    for (int k = 0; k < resolution; ++k) {
        for (int j = 0; j < resolution; ++j) {
            for (int i = 0; i < resolution; ++i) {
                // a. 将整数格点坐标映射到 [0, 1] 区间，再乘以频率
                float x = static_cast<float>(i) / resolution * frequency;
                float y = static_cast<float>(j) / resolution * frequency;
                float z = static_cast<float>(k) / resolution * frequency;

                // b. 获取该点的噪声值
                // Perlin噪声输出范围约[-1,1]，我们先把它映射到[0,1]
                float noise_value = perlin.getValue(x, y, z);
                float density = (noise_value + 1.0f) / 2.0f; 

                // c. 将密度值写入MarchingCubes对象
                mc.set_data(density, i, j, k);
                }
            }
        }
    //5. 运行MarchingCubes算法
    mc.run(threshold);
    //6. 将生成的网格写入PLY文件
    const char* filename = "cloud_mesh.ply";
    mc.writePLY(filename, false);

    std::cout << "over" << std::endl;
    // 释放资源
    mc.clean_all();

    return 0;
}