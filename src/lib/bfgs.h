/**
 * @file bfgs.h
 * @brief BFGS优化算法实现

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

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "matrix.h"

/// 浮点数三角矩阵，用于存储BFGS的Hessian近似矩阵
typedef triangular_matrix<fl> flmat;

/// @brief 计算负的矩阵向量乘积，用于BFGS算法中的搜索方向更新
template<typename Change>
void minus_mat_vec_product(const flmat& m, const Change& in, Change& out) {
	sz n = m.dim();
	VINA_FOR(i, n) {
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(m.index_permissive(i, j)) * in(j);
		out(i) = -sum;
	}
}

/// @brief 计算两个向量点乘
/// @return a·b
template<typename Change>
inline fl scalar_product(const Change& a, const Change& b, sz n) {
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a(i) * b(i);
	return tmp;
}

/// @brief Hessian矩阵更新
template<typename Change>
inline bool bfgs_update(flmat& h, const Change& p, const Change& y, const fl alpha) {
	const fl yp  = scalar_product(y, p, h.dim());
	if(alpha * yp < epsilon_fl) return false; // FIXME?
	Change minus_hy(y); minus_mat_vec_product(h, y, minus_hy);
	const fl yhy = - scalar_product(y, minus_hy, h.dim());
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	const sz n = p.num_floats();
	VINA_FOR(i, n)
		VINA_RANGE(j, i, n) // includes i
			h(i, j) +=   alpha * r * (minus_hy(i) * p(j)
	                                + minus_hy(j) * p(i)) +
			           + alpha * alpha * (r*r * yhy  + r) * p(i) * p(j); // s * s == alpha * alpha * p * p
	return true;
}

/// @brief 线性搜索
template<typename F, typename Conf, typename Change>
fl line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0, const Change& p, Conf& x_new, Change& g_new, fl& f1, int& evalcount) { // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	const fl pg = scalar_product(p, g, n);

	VINA_U_FOR(trial, max_trials) {
		// 计算f(x_k + alpha * p_k)
		x_new = x; x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		evalcount++;
		// 判断是否满足Armijo条件
		if(f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;	// 将步长减小
	}
	return alpha;
}

/// @brief 设置浮点数矩阵的对角线元素为指定值
inline void set_diagonal(flmat& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

/// @brief 从b中减去a，结果存储在b中，即计算两个变化量的差
template<typename Change>
void subtract_change(Change& b, const Change& a, sz n) { // b -= a
	VINA_FOR(i, n)
		b(i) -= a(i);
}

/// @brief BFGS拟牛顿优化算法主函数
/// @param f 目标函数对象，quasi_newton_aux::operator()
/// @param x 初始构象		初始位置（输入输出），优化后的最优位置
/// @param g 构象变化		梯度向量（输入输出）
/// @param max_steps 最大迭代步数
/// @param average_required_improvement 平均改进要求（未使用）
/// @param over 窗口大小（未使用）
/// @param evalcount 函数评估次数计数器
/// @return 最优函数值
template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const unsigned max_steps, const fl average_required_improvement, const sz over,
		int& evalcount) { // x is I/O, final value is returned
	sz n = g.num_floats();
	// 初始化Hessian矩阵的逆矩阵近似为单位矩阵
	flmat h(n, 0);
	set_diagonal(h, 1);

	// 复制初始梯度和构象
	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);
	evalcount++;

	// 原始状态存档
	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);

	Change p(g);	// 搜索方向，用梯度初始化

	// 结合能历史记录，似乎只记录但是没有使用
	flv f_values; f_values.reserve(max_steps+1);
	f_values.push_back(f0);

	VINA_U_FOR(step, max_steps) {
		// 计算搜索方向  p = -H * g
		minus_mat_vec_product(h, g, p);

		// 线性搜索确定最优步长，并更新构象x_new
		fl f1 = 0;
		const fl alpha = line_search(f, n, x, g, f0, p, x_new, g_new, f1, evalcount);

		// 计算梯度变化，y=g_new - g
		Change y(g_new); subtract_change(y, g, n);

		// 更新能量与构象
		f_values.push_back(f1);
		f0 = f1;
		x = x_new;

		// 检查收敛条件，若满足要求则说明到达一个梯度很小的平稳点（很可能就是局部极小值点），迭代终止
		if(!(std::sqrt(scalar_product(g, g, n)) >= 1e-5)) break; // breaks for nans too // FIXME !!?? 
		g = g_new; // ?

		// 启发式的初始化Hessian矩阵的策略，仅在第一次迭代时应用
		if(step == 0) {
			const fl yy = scalar_product(y, y, n);
			if(std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * scalar_product(y, p, n) / yy);
		}

		// 更新Hessian逆矩阵近似
		bool h_updated = bfgs_update(h, p, y, alpha);
	}
	// 如果优化失败（函数值增加或出现NaN），恢复到原始状态
	if(!(f0 <= f_orig)) { // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

#endif
