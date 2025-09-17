#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "NOISE/PerlinNoise.hpp"
#include "Surface/MarchingCubes.h"
#include "Grid_Construction/vector3D.hpp" // 确保包含了 Vector3D

// 平滑过渡函数
float smoothstep(float edge0, float edge1, float x) {
    x = std::max(0.0f, std::min(1.0f, (x - edge0) / (edge1 - edge0)));
    return x * x * (3.0f - 2.0f * x);
}

// 定义一个云团块的结构体
struct CloudClump {
    Vector3D center; // 云团中心位置
    float radius;    // 云团半径
};

int main() {
    // 1. 噪声生成器
    unsigned int seed = 1337;
    PerlinNoise perlin(seed);
    std::cout << "使用种子: " << seed << "\n";

    // 2. 参数定义
    const int resolution = 128; 
    const float threshold = 0.25f;

    // 3. MarchingCubes 初始化
    MarchingCubes mc(resolution, resolution, resolution);
    mc.init_all();

    std::cout << "正在生成模型文件，分辨率 " << resolution << "^3 ...\n";

    // ==================== 坐标系修正后的核心代码 ====================

    // A. 定义组成积云的多个“云团块” (现在中心点Z坐标代表高度)
    std::vector<CloudClump> clumps;
    // 一个较大的主基座云团，位置较低
    clumps.push_back({Vector3D(0.0f, 0.0f, -0.2f), 0.5f});
    // 几个堆叠在上面的、小一些的、位置随机的云团
    clumps.push_back({Vector3D(0.2f, 0.1f, 0.0f), 0.35f});
    clumps.push_back({Vector3D(-0.25f, -0.1f, 0.1f), 0.4f});
    clumps.push_back({Vector3D(-0.0f, 0.2f, 0.2f), 0.3f});
    clumps.push_back({Vector3D(0.1f, -0.15f, 0.25f), 0.3f});


    // B. 噪声参数
    float shape_freq = 1.5f;
    int detail_octaves = 5;
    float detail_freq = 8.0f;
    float detail_persistence = 0.5f;
    
    // C. 物理启发参数 (现在基于Z轴)
    float cloud_scale_xy = 1.0f;       // 云朵在水平方向(XY)的拉伸
    float cloud_scale_z = 1.5f;        // 云朵在垂直方向(Z)的拉伸
    float cloud_center_z = -0.2f;      // 云朵中心的垂直位置 [-1, 1]
    float condensation_height = -0.6f; // **物理上的“凝结高度” (Z轴)**
    float condensation_falloff = 0.25f;

    // 4. 核心循环
    for (int k = 0; k < resolution; ++k) {
        for (int j = 0; j < resolution; ++j) {
            for (int i = 0; i < resolution; ++i) {
                
                // a. 坐标转换到 [-1, 1] 空间
                float x = 2.0f * (static_cast<float>(i) / (resolution - 1) - 0.5f);
                float y = 2.0f * (static_cast<float>(j) / (resolution - 1) - 0.5f);
                float z = 2.0f * (static_cast<float>(k) / (resolution - 1) - 0.5f);
                
                // b. 计算多个云团聚合后的基础密度
                float max_clump_density = 0.0f;
                for (const auto& clump : clumps) {
                    // **核心修正：距离计算时，Z轴应用垂直拉伸**
                    float dx = (x - clump.center.x) / cloud_scale_xy;
                    float dy = (y - clump.center.y) / cloud_scale_xy;
                    float dz = (z - clump.center.z) / cloud_scale_z;
                    float dist_sq = dx*dx + dy*dy + dz*dz;
                    
                    // 用半径的平方来定义影响范围
                    float clump_radius_sq = (clump.radius * clump.radius);

                    if (dist_sq < 1.0f) {
                        float falloff = 1.0f - smoothstep(0.0f, 1.0f, dist_sq);
                        max_clump_density = std::max(max_clump_density, falloff);
                    }
                }
                
                // c. 用低频噪声侵蚀大的形状
                float shape_noise = (perlin.getValue(x * shape_freq, y * shape_freq, z * shape_freq) + 1.0f) / 2.0f;
                float eroded_density = max_clump_density - (1.0f - shape_noise);
                eroded_density = std::max(0.0f, eroded_density);

                // d. **核心修正：模拟物理平底 (作用于Z轴)**
                float bottom_mask = smoothstep(condensation_height, condensation_height + condensation_falloff, z);
                float final_density = eroded_density * bottom_mask;

                // e. 添加高频细节
                if (final_density > 0.01f) {
                    float detail_noise = 0.0f;
                    float freq = detail_freq;
                    float amp = 1.0f;
                    for (int l = 0; l < detail_octaves; ++l) {
                        detail_noise += perlin.getValue(x * freq, y * freq, z * freq) * amp;
                        freq *= 2.0f;
                        amp *= detail_persistence;
                    }
                    final_density -= detail_noise * detail_noise * 0.1f * smoothstep(0.0f, 0.3f, final_density);
                    final_density = std::max(0.0f, final_density);
                }

                // f. 写入数据
                mc.set_data(final_density, i, j, k);
            }
        }
    }

    // =================================================================

    // 5. 运行Marching Cubes
    mc.run(threshold);

    // 6. 保存结果
    const char* filename = "cumulus_cloud_v3_correct_axis.ply";
    mc.writePLY(filename, false);

    std::cout << "成功生成: " << filename << std::endl;
    
    // 7. 清理
    mc.clean_all();

    return 0;
}