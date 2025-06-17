/**
 * @file grid_dim.h
 * @brief 网格维度定义

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_GRID_DIM_H
#define VINA_GRID_DIM_H

#include <boost/array.hpp>

#include "common.h"

/**
 * @struct grid_dim
 * @brief 定义了一维网格的起始点、结束点和体素数量
 */
struct grid_dim {
    fl begin;      ///< 网格起始坐标
    fl end;        ///< 网格结束坐标
    sz n_voxels;   ///< 间隔区间数
    
    /**
     * @brief 默认构造函数初始化所有成员为0
     */
    grid_dim() : begin(0), end(0), n_voxels(0) {}
    
    /**
     * @brief 计算网格跨度
     * @return 网格的跨度长度（end - begin）
     */
    fl span() const { return end - begin; }
    
    /**
     * @brief 检查网格是否启用
     * @return 如果体素数量大于0则返回true，否则返回false
     */
    bool enabled() const { return (n_voxels > 0); }
    
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
        ar & begin;
        ar & end;
        ar & n_voxels;
    }
};

/**
 * @brief 比较两个grid_dim对象是否相等
 */
inline bool eq(const grid_dim& a, const grid_dim& b) {
    return a.n_voxels == b.n_voxels && eq(a.begin, b.begin) && eq(a.end, b.end);
}

/**
 * @typedef grid_dims
 * @brief 三维网格维度数组类型
 */
typedef boost::array<grid_dim, 3> grid_dims;

/**
 * @brief 比较两个grid_dims对象是否相等
 */
inline bool eq(const grid_dims& a, const grid_dims& b) {
    return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

/**
 * @brief 打印网格维度信息
 */
inline void print(const grid_dims& gd, std::ostream& out = std::cout) {
    VINA_FOR_IN(i, gd)  // 遍历三个维度
        std::cout << gd[i].n_voxels << " [" << gd[i].begin << " .. " << gd[i].end << "]\n";
}

/**
 * @brief 获取网格的起始坐标向量
 */
inline vec grid_dims_begin(const grid_dims& gd) {
    vec tmp;
    VINA_FOR_IN(i, gd)  // 遍历三个维度，提取begin值
        tmp[i] = gd[i].begin;
    return tmp;
}

/**
 * @brief 获取网格的结束坐标向量
 */
inline vec grid_dims_end(const grid_dims& gd) {
    vec tmp;
    VINA_FOR_IN(i, gd)  // 遍历三个维度，提取end值
        tmp[i] = gd[i].end;
    return tmp;
}

#endif
