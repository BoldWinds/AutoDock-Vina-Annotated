/**
 * @file coords.h
 * @brief 定义计算RMSD上界计算、查找最相似构象和添加对接结果到输出容器

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

#ifndef VINA_COORDS_H
#define VINA_COORDS_H

#include "conf.h"
#include "atom.h" // for atomv

/**
 * @brief 计算两组坐标向量间的RMSD上界
 */
fl rmsd_upper_bound(const vecv& a, const vecv& b);

/**
 * @brief 在输出容器中查找与给定坐标最相似的构象
 */
std::pair<sz, fl> find_closest(const vecv& a, const output_container& b);

/**
 * @brief 向输出容器添加新的对接结果
 * @param out 输出容器，存储对接结果
 * @param t 待添加的输出类型结果
 * @param min_rmsd 最小RMSD阈值，用于判断构象相似性
 * @param max_size 容器最大容量限制
 * @note 实现去重逻辑：相似构象保留能量更低者，容器满时替换最差结果
 */
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size);

#endif
