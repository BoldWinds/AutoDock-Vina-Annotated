/**
 * @file parallel_mc.h
 * @brief 并行蒙特卡洛搜索功能的头文件定义

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

#ifndef VINA_PARALLEL_MC_H
#define VINA_PARALLEL_MC_H

#include "monte_carlo.h"

/**
 * @brief 并行蒙特卡洛搜索类
 * @details 管理多个蒙特卡洛搜索任务的并行执行，用于提高分子对接的搜索效率
 */
struct parallel_mc {
	monte_carlo mc;
	sz num_tasks;
	sz num_threads;
	bool display_progress;

	parallel_mc() : num_tasks(8), num_threads(1), display_progress(true) {}

    /**
     * @brief 执行并行蒙特卡洛搜索
     * @param m 对接模型
     * @param out 输出结果容器
     * @param p 按原子预计算数据
     * @param ig 刚体能量网格
     * @param corner1 搜索空间第一个角点
     * @param corner2 搜索空间第二个角点
     * @param generator 随机数生成器
     * @param progress_callback 进度回调函数指针，目前无作用
     */
	void operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig, const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback) const;
};

#endif
