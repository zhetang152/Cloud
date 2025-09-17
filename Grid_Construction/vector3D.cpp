#include"Vector3D.hpp"
#include <cassert>

Vector3D Vector3D::operator+(const Vector3D& other)const{
    return Vector3D(x+other.x,y+other.y,z+other.z);
}
Vector3D Vector3D::operator-(const Vector3D& other)const{
    return Vector3D(x-other.x,y-other.y,z-other.z);
}
Vector3D Vector3D::operator*(float scalar)const{
    return Vector3D(x*scalar,y*scalar,z*scalar);
}
Vector3D Vector3D::operator/(float scalar) const{
    assert(std::abs(scalar) > 1e-9f);
    return (*this) * (1.0f / scalar);
}

float Vector3D::length() const{
    return std::sqrt(x*x+y*y+z*z);
}
Vector3D& Vector3D::normalize(){
    float len = this->length();
    if(len>1e-6f){
        this->x/=len;
        this->y/=len;
        this->z/=len;
    }
    return *this;
}

float Vector3D::dot(const Vector3D& other)const{
    return Vector3D::x*other.x+Vector3D::y*other.y+Vector3D::z*other.z;
}
Vector3D operator*(float scalar, const Vector3D& vec) {
    // 直接调用已有的成员函数，实现代码复用
    return vec * scalar;
}