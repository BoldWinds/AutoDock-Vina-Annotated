/**
 * @file szv_grid.h
 * @brief svz网格定义，主要功能优化网格维度和筛选可能相互作用的原子

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

#ifndef VINA_SZV_GRID_H
#define VINA_SZV_GRID_H

#include "model.h"
#include "grid_dim.h"
#include "array3d.h"

/**
 * @brief SZV(Size Vector)网格结构 - 空间分割索引系统
 * 
 * 将3D空间划分为网格体素，每个体素存储在截断距离内的原子索引列表，
 * 用于优化分子对接过程中的距离查询性能
 */
struct szv_grid {
    szv_grid() {}
    
    /**
     * @brief 构造函数，根据分子模型和网格维度初始化SZV网格
     */
    szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr);
    
    /**
     * @brief 根据坐标获取可能相互作用的原子索引列表
     * @param coords 查询的3D坐标点
     */
    const szv& possibilities(const vec& coords) const;
    
    /**
     * @brief 计算每个网格体素平均包含的原子数量
     */
    fl average_num_possibilities() const;

private:
    array3d<szv> m_data;  ///< 三维数组，每个元素存储原子索引向量(szv)
    vec m_init;				    ///< 网格起始坐标
    vec m_range;			    ///< 网格空间范围（各轴的长度）
    
    /**
     * @brief 将网格索引转换为实际坐标
     */
    vec index_to_coord(sz i, sz j, sz k) const;
};

/**
 * @brief 为SZV网格生成优化的网格维度配置
 * @param gd 原始网格维度
 * @return 优化后的网格维度，区间距大小约为3Å
 */
grid_dims szv_grid_dims(const grid_dims& gd);


#endif
