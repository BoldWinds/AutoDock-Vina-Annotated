/**
 * @file parallel_mc.cpp
 * @brief 并行蒙特卡洛搜索功能的实现

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

#include "parallel.h"
#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"

/**
 * @brief 并行蒙特卡洛任务结构
 * @details 包含单个蒙特卡洛搜索任务所需的对接模型、对接结果、随机数生成器
 */
struct parallel_mc_task {
	model m;
	output_container out;
	rng generator;
	parallel_mc_task(const model& m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

/**
 * @brief 并行蒙特卡洛辅助执行器
 * @details 包含执行蒙特卡洛搜索所需的所有共享参数和执行逻辑
 */
struct parallel_mc_aux {
	const monte_carlo* mc;
	const precalculate_byatom* p;
	const igrid* ig;
	const vec* corner1;
	const vec* corner2;
	parallel_progress* pg;
	parallel_mc_aux(const monte_carlo* mc_, const precalculate_byatom* p_, const igrid* ig_, const vec* corner1_, const vec* corner2_, parallel_progress* pg_)
		: mc(mc_), p(p_), ig(ig_), corner1(corner1_), corner2(corner2_), pg(pg_) {}
	// 执行monte-carlo搜索
	void operator()(parallel_mc_task& t) const {
		(*mc)(t.m, t.out, *p, *ig, *corner1, *corner2, pg, t.generator);
	}
};

/**
 * @brief 将in中的output_type放入out中
 */
void merge_output_containers(const output_container& in, output_container& out, fl min_rmsd, sz max_size) {
	VINA_FOR_IN(i, in)
		add_to_output_container(out, in[i], min_rmsd, max_size);
}

/**
 * @brief 合并多个并行任务的输出容器到out中，返回的out是按照结合能排序的
 */
void merge_output_containers(const parallel_mc_task_container& many, output_container& out, fl min_rmsd, sz max_size) {
	min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
	VINA_FOR_IN(i, many)
		merge_output_containers(many[i].out, out, min_rmsd, max_size);
	out.sort();
}

/**
 * @brief 执行并行蒙特卡洛搜索的主要函数
 * @param m 分子模型
 * @param out 输出结果容器
 * @param p 原子间预计算数据
 * @param ig 相互作用网格
 * @param corner1 搜索空间第一个角点
 * @param corner2 搜索空间第二个角点
 * @param generator 随机数生成器
 * @param progress_callback 进度回调函数指针
 * @note 创建多个独立的蒙特卡洛任务，使用不同的随机种子并行执行，最后合并所有结果
 */
void parallel_mc::operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig, const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback) const {
	// 初始化进度管理器
	parallel_progress pp (progress_callback);
	// 创建辅助执行器，传递所有必要的参数
	parallel_mc_aux parallel_mc_aux_instance(&mc, &p, &ig, &corner1, &corner2, (display_progress ? (&pp) : NULL));
	// 创建任务容器并填充任务
	parallel_mc_task_container task_container;
	VINA_FOR(i, num_tasks)
		task_container.push_back(new parallel_mc_task(m, random_int(0, 1000000, generator)));
	// 如果需要显示进度，初始化进度计数器
	if(display_progress) 
		pp.init(num_tasks * mc.global_steps);
	// 创建并行迭代器并执行所有任务
	parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task, true> parallel_iter_instance(&parallel_mc_aux_instance, num_threads);
	parallel_iter_instance.run(task_container);
	// 合并所有任务的结果到最终输出容器
	merge_output_containers(task_container, out, mc.min_rmsd, mc.num_saved_mins);
}
