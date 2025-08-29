#pragma once
#include<iostream>
#include<vector>
#include<numeric>
#include<algorithm>
#include<chrono>
#include<string>
#include<random>
#include<limits>

/**
 * @brief 3D整数坐标的哈希函数
 * @param x,y,z 整数坐标
 * @param seed 用于初始化Hash过程的种子, 允许生成不同的Hash结果
 * @return 返回一个确定性的伪随机无符号整数
*/
inline unsigned int hash(int x, int y, int z, unsigned int seed) {
    //将有符号整数转换为无符号整数, 以便进行位运算
    unsigned int ux = static_cast<unsigned int>(x);
    unsigned int uy = static_cast<unsigned int>(y);
    unsigned int uz = static_cast<unsigned int>(z);
    
    // 使用传入的种子作为哈希状态的初始值
    unsigned int h = seed;

    //选择一些大的素数作为乘子
    //以便最大化雪崩效应和分布均匀性
    const unsigned int prime1 = 73856093;
    const unsigned int prime2 = 19349663;
    const unsigned int prime3 = 83492791;
    
    //混合过程
    h = (h ^ (ux * prime1));
    h = (h ^ (uy * prime2));
    h = (h ^ (uz * prime3));
    
    //最终混合
    h = h ^ (h >> 16);
    h = h * 0x85ebca6b;
    h = h ^ (h >> 13);
    h = h * 0xc2b2ae35;
    h = h ^ (h >> 16);

    return h;
}