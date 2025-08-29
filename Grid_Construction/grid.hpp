#pragma once

#include<vector>
#include<cassert>
#include<algorithm>
template<typename T>
class Grid{
private:
    int m_width;
    int m_height;
    int m_depth;
    std::vector<T> m_data;
public:
    Grid(int width, int height,int depth,const T& initial_value=T()):m_width(width), m_height(height), m_depth(depth){
        m_data.resize(width*height*depth, initial_value);
        }
    T& operator()(int i, int j, int k) {
        // 使用三维索引公式
        assert(i >= 0 && i < m_width  && "Grid3D access out of bounds (i)");
        assert(j >= 0 && j < m_height && "Grid3D access out of bounds (j)");
        assert(k >= 0 && k < m_depth  && "Grid3D access out of bounds (k)");
        return m_data[i + j * m_width + k * m_width * m_height];
    }
    const T& operator()(int i, int j, int k) const {
        assert(i >= 0 && i < m_width  && "Grid3D access out of bounds (i)");
        assert(j >= 0 && j < m_height && "Grid3D access out of bounds (j)");
        assert(k >= 0 && k < m_depth  && "Grid3D access out of bounds (k)");
        return m_data[i + j * m_width + k * m_width * m_height];
    }
    void fill(const T& value) {
        std::fill(m_data.begin(), m_data.end(), value);
    }
    int getWidth() const { return m_width; }
    int getHeight() const { return m_height; }
    int getDepth() const { return m_depth; }
};