/**
 * @file random.h
 * @brief 随机数生成相关函数声明

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

#ifndef VINA_RANDOM_H
#define VINA_RANDOM_H

#include <random>
#include <boost/random.hpp>
#include "common.h"

/**
 * @brief 随机数生成器类型定义，使用Boost的梅森旋转算法
 */
typedef boost::mt19937 rng;

/**
 * @brief 生成指定范围内的随机浮点数
 * @param a 范围下界（不包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回(a, b]范围内的随机浮点数
 * @pre 要求 a < b
 */
fl random_fl(fl a, fl b, rng& generator);

/**
 * @brief 生成符合正态分布的随机浮点数
 * @param mean 正态分布的均值
 * @param sigma 正态分布的标准差
 * @param generator 随机数生成器引用
 * @return 返回符合正态分布的随机浮点数
 * @pre 要求 sigma >= 0
 */
fl random_normal(fl mean, fl sigma, rng& generator);

/**
 * @brief 生成指定范围内的随机整数
 * @param a 范围下界（包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回[a, b]范围内的随机整数
 * @pre 要求 a <= b
 */
int random_int(int a, int b, rng& generator);

/**
 * @brief 生成指定范围内的随机无符号整数
 * @param a 范围下界（包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回[a, b]范围内的随机无符号整数
 * @pre 要求 a <= b
 */
sz random_sz(sz a, sz b, rng& generator);

/**
 * @brief 生成单位球体内的随机三维向量
 * @param generator 随机数生成器引用
 * @return 返回以原点为中心、半径为1的球体内的随机三维向量
 * @note 使用拒绝采样方法，平均需要运行约2次
 */
vec random_inside_sphere(rng& generator);

/**
 * @brief 生成指定矩形框内的随机三维向量
 * @param corner1 矩形框的第一个角点坐标
 * @param corner2 矩形框的第二个角点坐标
 * @param generator 随机数生成器引用
 * @return 返回矩形框内的随机三维向量
 * @pre 要求 corner1[i] < corner2[i] 对所有维度i成立
 */
vec random_in_box(const vec& corner1, const vec& corner2, rng& generator);

/**
 * @brief 自动生成随机种子
 * @return 返回基于系统熵和时间的随机种子值
 * @note 使用std::random_device获取系统熵，结合64位梅森旋转算法生成种子
 */
int auto_seed();

#endif
