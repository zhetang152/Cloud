#pragma once
#include"vector3D.hpp"
#include"grid.hpp"

enum class CellType{
    AIR,
    FLUID,
    SOLID
};
struct Particles{
    /* data */
    Vector3D position;
    Vector3D velocity;
};
class MACGrid{
private:
    int x,y,z;
    float dx;
    Grid<CellType> m_celltypes;
    Grid<float> m_u;
    Grid<float> m_v;
    Grid<float> m_w;
    Grid<float> m_pressure;
    Grid<Vector3D> m_solidvelocity;
    Grid<float> m_liquid_phi; //液体的SDF
    Grid<float> m_volume_fractions;//体积分数
    Grid<float> m_area_u, m_area_v, m_area_w; //用于计算速度的开放面积分数
    std::vector<Particles> m_particles; //用于标记粒子位置的向量
    Grid<float> m_density;
    Grid<float> m_temperature;
    public:
    MACGrid(int tx, int ty, int tz, float initial_dx):
    x(tx),y(ty),z(tz),dx(initial_dx),
    m_celltypes(tx,ty,tz,CellType::AIR),
    m_pressure(tx,ty,tz),
    m_solidvelocity(tx, ty, tz, Vector3D(0,0,0)),
    m_u(tx + 1, ty, tz),
    m_v(tx, ty + 1, tz),
    m_w(tx, ty, tz + 1),
    m_liquid_phi(tx, ty, tz, 0.0f),
    m_volume_fractions(tx, ty, tz, 0.0f),
    m_area_u(tx + 1, ty, tz, 0.0f),
    m_area_v(tx, ty + 1, tz, 0.0f),
    m_area_w(tx, ty, tz + 1, 0.0f),
    m_density(tx,ty,tz,0.0f),
    m_temperature(tx,ty,tz,0.0f)
    {
    }
    //提供访问
    int getDimX() const { return x; };
    int getDimY() const { return y; };
    int getDimZ() const { return z; };
    Grid<float>& pressure() {return m_pressure;}
    //使用例: MAC.pressure()(i, j, k) = 10.f;
    const Grid<float>& pressure() const {return m_pressure;}
    
    Grid<Vector3D>& solidvelocity() { return m_solidvelocity; }
    const Grid<Vector3D>& solidvelocity() const { return m_solidvelocity; }
    std::vector<Particles>& particles() { return m_particles; }
    const std::vector<Particles>& particles() const { return m_particles; }
    Grid<float>& liquid_phi() {return m_liquid_phi;}
    const Grid<float>& liquid_phi() const {return m_liquid_phi;}
    Grid<float>& volumeFractions() { return m_volume_fractions; }
    const Grid<float>& volumeFractions() const { return m_volume_fractions; }
    Grid<float>& area_u() { return m_area_u; }
    const Grid<float>& area_u() const { return m_area_u; }
    Grid<float>& area_v() { return m_area_v; }
    const Grid<float>& area_v() const { return m_area_v; }
    Grid<float>& area_w() { return m_area_w; }
    const Grid<float>& area_w() const { return m_area_w; }
    Grid<float>& u() { return m_u; }
    const Grid<float>& u() const { return m_u; }
    Grid<float>& u() { return m_u; }
    const Grid<float>& u() const { return m_u; }
    Grid<float>& v() { return m_v; }
    const Grid<float>& v() const { return m_v; }
    Grid<float>& w() { return m_w; }
    const Grid<float>& w() const { return m_w; }
    Grid<float>& density() { return m_density; }
    const Grid<float>& density() const { return m_density; }
    Grid<float>& temperature() { return m_temperature; }
    const Grid<float>& temperature() const { return m_temperature; }
    Grid<CellType>& celltypes() { return m_celltypes; }
    const Grid<CellType>& celltypes() const { return m_celltypes; }
    float getDx()const{return dx;}
    Vector3D PositionOfPressure(int i, int j,int k) const{
        return Vector3D(
            (static_cast<float>(i)+0.5f)*dx,
            (static_cast<float>(j) + 0.5f) * dx,
            (static_cast<float>(k) + 0.5f) * dx
        );
    }
    Vector3D positionOfU(int i, int j, int k) const {
        return Vector3D(
            static_cast<float>(i) * dx,
            (static_cast<float>(j) + 0.5f) * dx,
            (static_cast<float>(k) + 0.5f) * dx
        );
    }
    Vector3D positionOfV(int i, int j, int k) const{
        return Vector3D(
            (static_cast<float>(i) + 0.5f) * dx,
            static_cast<float>(j) * dx,
            (static_cast<float>(k) + 0.5f) * dx
        );
    }
    Vector3D positionOfW(int i, int j, int k) const {
        return Vector3D(
            (static_cast<float>(i) + 0.5f) * dx,
            (static_cast<float>(j) + 0.5f) * dx,
            static_cast<float>(k) * dx
        );
    }
};