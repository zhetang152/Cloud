#include "GravityForce.hpp"

GravityForce::GravityForce(const Vector3D& gravity): m_gravity(gravity) {}
void GravityForce::apply(MACGrid& grid, float dt) {
    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j< grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                grid.v()(i, j, k) += m_gravity.y * dt;
            }
        }
    }
}