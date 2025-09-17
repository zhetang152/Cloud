#pragma once
#include<cmath>
#include<iostream>
#include "NoiseUtils.hpp"
#include "Grid_Construction\vector3D.hpp"

//先写一个线性插值
inline float lerp(float a, float b, float t){
    return a + t * (b - a);
}

class PerlinNoise{
private:
    unsigned int m_seed;
    static const Vector3D grad_table[12];
    //缓和曲线 6t^5 - 15t^4 + 10t^3
    float fade(float t) const{
        return t * t * t * (t * (t * 6 - 15) + 10);
    }
public:
    PerlinNoise(unsigned int seed):m_seed(seed){}
    //计算噪声值的主函数
    float getValue(float x, float y, float z) const{
        //向下取整, 得到基点坐标
        int X = static_cast<int>(std::floor(x));
        int Y = static_cast<int>(std::floor(y));
        int Z = static_cast<int>(std::floor(z));
        //计算相对坐标
        float xf = x - X;
        float yf = y - Y;
        float zf = z - Z;
        //梯度向量
        //定义整数坐标
        int X1 = X + 1;
        int Y1 = Y + 1;
        int Z1 = Z + 1;
        //A. 对8个顶点坐标进行Hash运算
        unsigned int h000 = hash(X, Y, Z, m_seed);
        unsigned int h100 = hash(X1, Y, Z, m_seed);
        unsigned int h010 = hash(X, Y1, Z, m_seed);
        unsigned int h110 = hash(X1, Y1, Z, m_seed);
        unsigned int h001 = hash(X, Y, Z1, m_seed);
        unsigned int h101 = hash(X1, Y, Z1, m_seed);
        unsigned int h011 = hash(X, Y1, Z1, m_seed);
        unsigned int h111 = hash(X1, Y1, Z1, m_seed);
        //B. 计算梯度向量
        const Vector3D& g000 = grad_table[h000 % 12];
        const Vector3D& g100 = grad_table[h100 % 12];
        const Vector3D& g010 = grad_table[h010 % 12];
        const Vector3D& g110 = grad_table[h110 % 12];
        const Vector3D& g001 = grad_table[h001 % 12];
        const Vector3D& g101 = grad_table[h101 % 12];
        const Vector3D& g011 = grad_table[h011 % 12];
        const Vector3D& g111 = grad_table[h111 % 12];
        //计算各顶点到采样点的影响
        //A. 计算各顶点到采样点的距离向量
        Vector3D d000(xf, yf, zf);
        Vector3D d100(xf - 1, yf, zf);
        Vector3D d010(xf, yf - 1, zf);
        Vector3D d110(xf - 1, yf - 1, zf);
        Vector3D d001(xf, yf, zf - 1);
        Vector3D d101(xf - 1, yf, zf - 1);
        Vector3D d011(xf, yf - 1, zf - 1);
        Vector3D d111(xf - 1, yf - 1, zf - 1);
        //计算点积
        float n000 = g000.dot(d000);
        float n100 = g100.dot(d100);
        float n010 = g010.dot(d010);
        float n110 = g110.dot(d110);
        float n001 = g001.dot(d001);
        float n101 = g101.dot(d101);
        float n011 = g011.dot(d011);
        float n111 = g111.dot(d111);
        //平滑插值
        //A. 计算缓和曲线值
        float u = fade(xf);
        float v = fade(yf);
        float w = fade(zf);
        //B. 三线性插值
        float ix0 = lerp(n000, n100, u);
        float ix1 = lerp(n010, n110, u);
        float ix2 = lerp(n001, n101, u);
        float ix3 = lerp(n011, n111, u);
        float iy0 = lerp(ix0, ix1, v);
        float iy1 = lerp(ix2, ix3, v);
        float value = lerp(iy0, iy1, w);
        return value;
    }
};

const Vector3D PerlinNoise::grad_table[12] = {
    Vector3D(1,1,0), Vector3D(-1,1,0), Vector3D(1,-1,0), Vector3D(-1,-1,0),
    Vector3D(1,0,1), Vector3D(-1,0,1), Vector3D(1,0,-1), Vector3D(-1,0,-1),
    Vector3D(0,1,1), Vector3D(0,-1,1), Vector3D(0,1,-1), Vector3D(0,-1,-1)
};