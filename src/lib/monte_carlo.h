/**
 * @file monte_carlo.h
 * @brief 蒙特卡洛算法定义

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

#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "incrementable.h"
#include "model.h"

struct monte_carlo {
    unsigned max_evals;         ///< 最大能量评估次数限制，0表示无限制
    unsigned global_steps;      ///< 全局搜索步数，控制蒙特卡罗采样次数
    fl temperature;             ///< 模拟温度，控制Metropolis接受概率
    vec hunt_cap;               ///< 局部优化的搜索范围限制向量
    fl min_rmsd;                ///< 最小可接受的RMSD值，用于去除相似构象
    sz num_saved_mins;          ///< 保存的最低能量构象数量
    fl mutation_amplitude;      ///< 构象变异的幅度参数
    unsigned local_steps;       ///< 局部优化步数，控制准牛顿法迭代次数

    // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  global_steps = 50*lig_atoms = 2500
    /**
     * @brief 默认构造函数
     * @note T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2
     * @note 除了mutation_amplitude和temperature，在初始化parallem_mc时都会根据具体参数初始化
     */
    monte_carlo() : max_evals(0), global_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) {}

	/// @brief  执行蒙特卡罗优化，返回单个最优结果
	output_type operator()(model& m, const precalculate_byatom& p, const igrid& ig, const vec& corner1,
                           const vec& corner2, incrementable* increment_me, rng& generator) const;
	/// @brief  执行蒙特卡罗优化，返回多个排序后的优化结果
	void operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                    const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
};

#endif

