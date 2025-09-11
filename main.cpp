#include<algorithm>
#include<chrono>
#include <fstream> // 用于文件读写
#include<iostream>
#include<vector>
#include "NOISE\PerlinNoise.hpp"
#include "Surface\MarchingCubes.h"

float smoothstep(float edge0, float edge1, float x){
    x = std::max(0.0f,std::min(1.0f, (x - edge0)/(edge1 - edge0)));
    return x * x * (3.0f - 2.0f * x);
}

int main() {
    // 1. 使用种子创建Perlin噪声生成器
    unsigned int seed = 1338; // 使用一个固定种子，方便每次生成同样的结果
    // unsigned int seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    PerlinNoise perlin(seed);
    std::cout << "使用种子: " << seed << "\n";

    // 2. 定义生成的参数
    const int resolution = 64; // 在一个 64x64x64 的立方体中采样
    //const float frequency = 5.0f;    // 噪声频率，控制云朵特征的大小
    const float threshold = 0.5f;   // 密度阈值, 作为marching cubes的iso值

    //3. 初始化MarchingCubes对象
    MarchingCubes mc(resolution, resolution, resolution);
    mc.init_all();

    std::cout << "正在生成OBJ文件，分辨率 " << resolution << "^3 ...\n";
    // fbm 参数
    int octaves = 6; //叠加层数
    float frequency = 2.0f; //初始频率
    float amplitude = 1.0f; //初始振幅
    float lacunarity = 2.0f; //频率倍增因子
    float persistence = 0.5f; //振幅衰减因子

    //云层形状参数
    float cloud_center_y = 0.6f; //云层中心高度
    float height = 0.4f; //云层垂直半径
    float bottom_falloff_start = 0.4f; //云层底部开始衰减的高度
    // 4. 遍历采样区域的每一个点
    for (int k = 0; k < resolution; ++k) {
        for (int j = 0; j < resolution; ++j) {
            for (int i = 0; i < resolution; ++i) {
                // a. 将整数格点坐标映射到 [0, 1] 区间，再乘以频率
                float x = static_cast<float>(i) / resolution;
                float y = static_cast<float>(j) / resolution;
                float z = static_cast<float>(k) / resolution;

                // b. 获取该点的fbm噪声值
                float total_noise = 0.0f;
                float current_freq = frequency;
                float current_amp = amplitude;
                float max_amp = 0.0f; // 用于归一化

                for(int l = 0; l<octaves; ++l){
                    total_noise += perlin.getValue(x * current_freq, y * current_freq, z * current_freq) * current_amp;
                    max_amp += current_amp;
                    current_freq *= lacunarity;
                    current_amp *= persistence;
                }
                float noise_value = total_noise / max_amp;

                //Perlin噪声输出范围约[-1,1]，我们先把它映射到[0,1]
                float density = (noise_value + 1.0f) / 2.0f; 
                // c.1 模拟热气泡的宏观形状
                float dist_from_center_sq = x*x + (y-0.2f)*(y-0.2f)*2.0f + z*z; // y方向压扁，中心稍微抬高
                float shape_mask = 1.0f - smoothstep(0.5f, 1.0f, dist_from_center_sq);
                // c.2 创建平坦云底
                float bottom_mask = smoothstep(-0.8f, -0.6f, y);

                // d. 最终密度值 = 细节噪声 * 形状遮罩
                float final_density = density * shape_mask * bottom_mask;

                // e. 将密度值写入MarchingCubes对象
                mc.set_data(density, i, j, k);
                }
            }
        }
    //5. 运行MarchingCubes算法
    mc.run(threshold);
    //6. 将生成的网格写入PLY文件
    const char* filename = "D:\\Code\\XLAB\\Fluid\\result\\cloud_mesh.ply";
    mc.writePLY(filename, false);

    std::cout << "over" << std::endl;
    // 释放资源
    mc.clean_all();

    return 0;
}