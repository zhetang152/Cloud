#include "advection.hpp"
#include <algorithm>

namespace Advector {

    // --- 底层插值辅助函数 ---

    float catmull_rom_1d(float p0, float p1, float p2, float p3, float t) {
        float t2 = t * t;
        float t3 = t2 * t;
        return p0 * (-0.5f * t + t2 - 0.5f * t3) +
               p1 * (1.0f - 2.5f * t2 + 1.5f * t3) +
               p2 * (0.5f * t + 2.0f * t2 - 1.5f * t3) +
               p3 * (-0.5f * t2 + 0.5f * t3);
    }

    float sample_linear(const Grid<float>& grid, const Vector3D& pos, float dx) {
        Vector3D grid_pos = pos / dx;

        int i = static_cast<int>(grid_pos.x);
        int j = static_cast<int>(grid_pos.y);
        int k = static_cast<int>(grid_pos.z);
        
        float tx = grid_pos.x - i;
        float ty = grid_pos.y - j;
        float tz = grid_pos.z - k;

        i = std::clamp(i, 0, grid.getWidth() - 2);
        j = std::clamp(j, 0, grid.getHeight() - 2);
        k = std::clamp(k, 0, grid.getDepth() - 2);

        float w000 = (1 - tx) * (1 - ty) * (1 - tz);
        float w100 = tx * (1 - ty) * (1 - tz);
        float w010 = (1 - tx) * ty * (1 - tz);
        float w110 = tx * ty * (1 - tz);
        float w001 = (1 - tx) * (1 - ty) * tz;
        float w101 = tx * (1 - ty) * tz;
        float w011 = (1 - tx) * ty * tz;
        float w111 = tx * ty * tz;

        return w000 * grid(i, j, k) + w100 * grid(i + 1, j, k) +
               w010 * grid(i, j + 1, k) + w110 * grid(i + 1, j + 1, k) +
               w001 * grid(i, j, k + 1) + w101 * grid(i + 1, j, k + 1) +
               w011 * grid(i, j + 1, k + 1) + w111 * grid(i + 1, j + 1, k + 1);
    }

    float sample_tricubic(const Grid<float>& grid, const Vector3D& pos, float dx) {
        Vector3D grid_pos = pos / dx;
        
        float x = std::clamp(grid_pos.x, 1.0f, static_cast<float>(grid.getWidth()) - 2.0001f);
        float y = std::clamp(grid_pos.y, 1.0f, static_cast<float>(grid.getHeight()) - 2.0001f);
        float z = std::clamp(grid_pos.z, 1.0f, static_cast<float>(grid.getDepth()) - 2.0001f);

        int ix = static_cast<int>(x);
        int iy = static_cast<int>(y);
        int iz = static_cast<int>(z);

        float tx = x - ix;
        float ty = y - iy;
        float tz = z - iz;

        float x_interp[4][4];
        for (int k_offset = 0; k_offset < 4; ++k_offset) {
            for (int j_offset = 0; j_offset < 4; ++j_offset) {
                x_interp[k_offset][j_offset] = catmull_rom_1d(
                    grid(ix - 1, iy + j_offset - 1, iz + k_offset - 1),
                    grid(ix    , iy + j_offset - 1, iz + k_offset - 1),
                    grid(ix + 1, iy + j_offset - 1, iz + k_offset - 1),
                    grid(ix + 2, iy + j_offset - 1, iz + k_offset - 1),
                    tx);
            }
        }

        float y_interp[4];
        for (int k_offset = 0; k_offset < 4; ++k_offset) {
            y_interp[k_offset] = catmull_rom_1d(
                x_interp[k_offset][0], x_interp[k_offset][1], x_interp[k_offset][2], x_interp[k_offset][3], ty);
        }

        return catmull_rom_1d(y_interp[0], y_interp[1], y_interp[2], y_interp[3], tz);
    }

    Vector3D get_velocity_at(const MACGrid& grid, const Vector3D& pos) {
        float dx = grid.getDx();
        
        float u_val = sample_linear(grid.u(), Vector3D(pos.x, pos.y - 0.5f * dx, pos.z - 0.5f * dx), dx);
        float v_val = sample_linear(grid.v(), Vector3D(pos.x - 0.5f * dx, pos.y, pos.z - 0.5f * dx), dx);
        float w_val = sample_linear(grid.w(), Vector3D(pos.x - 0.5f * dx, pos.y - 0.5f * dx, pos.z), dx);

        return Vector3D(u_val, v_val, w_val);
    }

    // --- 公共接口函数 ---

    Grid<float> advect(const Grid<float>& q_old, const MACGrid& velocityGrid, float dt) {
        int nx = q_old.getWidth();
        int ny = q_old.getHeight();
        int nz = q_old.getDepth();
        Grid<float> q_new(nx, ny, nz);
        float dx = velocityGrid.getDx();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    Vector3D pos = Vector3D((i + 0.5f) * dx, (j + 0.5f) * dx, (k + 0.5f) * dx);
                    Vector3D velocity = get_velocity_at(velocityGrid, pos);
                    Vector3D departure_pos = pos - velocity * dt;
                    q_new(i, j, k) = sample_tricubic(q_old, departure_pos, dx);
                }
            }
        }
        return q_new;
    }

    void advect_velocity(MACGrid& grid, float dt) {
        Grid<float> u_old = grid.u();
        Grid<float> v_old = grid.v();
        Grid<float> w_old = grid.w();
        float dx = grid.getDx();

        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX() + 1; ++i) {
                    Vector3D u_pos = grid.positionOfU(i, j, k);
                    Vector3D vel = get_velocity_at(grid, u_pos);
                    Vector3D departure_pos = u_pos - vel * dt;
                    grid.u()(i, j, k) = sample_linear(u_old, departure_pos, dx);
                }
            }
        }
        
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY() + 1; ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    Vector3D v_pos = grid.positionOfV(i, j, k);
                    Vector3D vel = get_velocity_at(grid, v_pos);
                    Vector3D departure_pos = v_pos - vel * dt;
                    grid.v()(i, j, k) = sample_linear(v_old, departure_pos, dx);
                }
            }
        }

        for (int k = 0; k < grid.getDimZ() + 1; ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    Vector3D w_pos = grid.positionOfW(i, j, k);
                    Vector3D vel = get_velocity_at(grid, w_pos);
                    Vector3D departure_pos = w_pos - vel * dt;
                    grid.w()(i, j, k) = sample_linear(w_old, departure_pos, dx);
                }
            }
        }
    }
}