#pragma once
#include<cmath>

class Vector3D{
public:
    float x,y,z;
    //默认构造函数
    Vector3D(): x(0.0f),y(0.0f),z(0.0f){}
    //含参构造函数
    Vector3D(float x,float y,float z):x(x),y(y),z(z){}
    //运算符重载
    Vector3D operator+(const Vector3D& other) const;
    Vector3D operator-(const Vector3D& other) const;
    Vector3D operator*(float scalar) const;
    Vector3D operator/(float scalar) const;
    //成员函数
    //点积
    float dot(const Vector3D& other) const;
    //length()函数读取x,y,z计算长度
    float length() const;
    //单位化向量
    Vector3D& normalize();
    //  
};
Vector3D operator*(float scalar, const Vector3D& vec);