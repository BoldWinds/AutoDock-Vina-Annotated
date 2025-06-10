/**
 * @file precalculate.h
 * @brief AutoDock Vina中的预计算模块，用于优化分子对接过程中的能量评估

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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include "scoring_function.h"
#include "matrix.h"

//Forward declaration
struct model;

/**
 * @class precalculate_element
 * @brief 预计算元素类，存储特定原子类型对的能量查找表
 * 
 * 该类为每对原子类型预计算并存储能量值，提供快速查找和插值功能
 */
struct precalculate_element
{
private:
    flv fast; ///< 快速查找表，牺牲精度换取速度
    fl factor; ///< 距离平方到数组索引的转换因子
public:
    prv smooth; ///< 存储(能量, 导数)对的数组，用于精确插值
    /**
     * @brief 构造函数
     * @param n 查找表大小
     * @param factor_ 把距离平方转换为对应的查找表索引的转换因子
     */
    precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) { }
    
    /**
     * @brief 快速能量评估（不计算导数）
     * @param r2 距离的平方
     * @return 相互作用能量
     * @note 使用预计算的查找表进行O(1)时间复杂度的能量查询
     */
    fl eval_fast(fl r2) const{
        assert(r2 * factor < fast.size());
        sz i = sz(factor * r2); // r2预期小于cutoff_sqr，且cutoff_sqr * factor + 1 < n，因此不会溢出
        assert(i < fast.size());
        return fast[i];
    };
    
    /**
     * @brief 评估能量及其导数
     * @param r2 距离的平方
     * @return pair<能量值, 对距离r的导数>
     * @note 先把r2变换为索引i1。实际上r2的位置在i1和i2之间，这时用p1和p2插值
     */
    pr eval_deriv(fl r2) const{
        fl r2_factored = factor * r2;
        assert(r2_factored + 1 < smooth.size());
        sz i1 = sz(r2_factored);
        sz i2 = i1 + 1; // r2预期小于cutoff_sqr，且cutoff_sqr * factor + 1 < n，因此不会溢出
        assert(i1 < smooth.size());
        assert(i2 < smooth.size());
        fl rem = r2_factored - i1; // 计算插值的余数部分
        assert(rem >= -epsilon_fl);
        assert(rem < 1 + epsilon_fl);
        const pr &p1 = smooth[i1];
        const pr &p2 = smooth[i2];
        // 线性插值计算能量和导数
        fl e = p1.first + rem * (p2.first - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
    };
    
    /**
     * @brief 从smooth数组的第一个元素初始化其他数据
     * @param rs 距离数组
     * @note 计算导数并填充fast查找表，使用数值微分计算导数
     */
    void init_from_smooth_fst(const flv &rs){
        sz n = smooth.size();
        VINA_CHECK(rs.size() == n);
        VINA_CHECK(fast.size() == n);
        VINA_FOR(i, n) {
            // 计算导数dor
            fl &dor = smooth[i].second;
            if (i == 0 || i == n - 1)
                dor = 0; // 边界点导数设为0
            else {
                fl delta = rs[i + 1] - rs[i - 1];
                fl r = rs[i];
                // 使用中心差分法计算数值导数
                dor = (smooth[i + 1].first - smooth[i - 1].first) / (delta * r);
            }
            // 从smooth.first计算fast查找表，使用相邻点的平均值
            fl f1 = smooth[i].first;
            fl f2 = (i + 1 >= n) ? 0 : smooth[i + 1].first;
            fast[i] = (f2 + f1) / 2;
        }
    };
    
    /**
     * @brief 找到smooth数组中能量值最小的索引
     * @return 最小能量值对应的索引
     */
    sz min_smooth_fst() const{
        sz tmp = 0; // 如果smooth为空则返回0
        VINA_FOR_IN(i_inv, smooth) {
            sz i = smooth.size() - i_inv - 1; // 反向遍历
            if (i_inv == 0 || smooth[i].first < smooth[tmp].first)
                tmp = i;
        }
        return tmp;
    };
    
    /**
     * @brief 扩展smooth数组的第一个元素（能量值）
     * @param rs 距离数组
     * @param left 向左扩展的距离
     * @param right 向右扩展的距离
     * @note 在最优距离附近创建平坦区域，用于能量函数的平滑化
     */
    void widen_smooth_fst(const flv &rs, fl left, fl right){
        flv tmp(smooth.size(), 0); // 新的smooth[].first值
        sz min_index = min_smooth_fst();
        VINA_CHECK(min_index < rs.size()); // 对于n == 0不成立
        VINA_CHECK(rs.size() == smooth.size());
        fl optimal_r = rs[min_index]; // 最优距离
        VINA_FOR_IN(i, smooth)
        {
            fl r = rs[i];
            // 在最优距离附近创建平坦区域
            if (r < optimal_r - left)
                r += left;
            else if (r > optimal_r + right)
                r -= right;
            else
                r = optimal_r; // 在平坦区域内

            // 边界检查
            if (r < 0)
                r = 0;
            if (r > rs.back())
                r = rs.back();

            tmp[i] = eval_deriv(sqr(r)).first;
        }
        VINA_FOR_IN(i, smooth)
        smooth[i].first = tmp[i];
    };
    
    /**
     * @brief 扩展能量函数并重新初始化
     * @param rs 距离数组
     * @param left 向左扩展的距离
     * @param right 向右扩展的距离
     */
    void widen(const flv &rs, fl left, fl right){
        widen_smooth_fst(rs, left, right);
        init_from_smooth_fst(rs);
    };
};

/**
 * @class precalculate
 * @brief 基于原子类型的预计算类
 * 
 * 为所有原子类型对预计算相互作用能量，建立查找表用于快速能量评估
 */
struct precalculate 
{
private:
    /**
     * @brief 计算距离数组
     * @return 距离值数组
     * @note 根据索引计算对应的距离值：r = sqrt(i / factor)
     */
    flv calculate_rs() const{
        flv tmp(m_n, 0);
        VINA_FOR(i, m_n)
            tmp[i] = std::sqrt(i / m_factor);
        return tmp;
    };

    fl m_cutoff_sqr; ///< 截断距离的平方
    fl m_max_cutoff_sqr; ///< 最大截断距离的平方
    sz m_n; ///< 查找表大小
    fl m_factor; ///< 距离平方到索引的转换因子

    triangular_matrix<precalculate_element> m_data; ///< 存储所有原子类型对预计算数据的三角矩阵
public:
    /**
     * @brief 默认构造函数
     */
    precalculate() { }
    
    /**
     * @brief 构造函数，基于评分函数初始化预计算表
     * @param sf 评分函数对象
     * @param v 能量上限值，默认为最大浮点数
     * @param factor 距离平方到索引的转换因子，默认为32
     * @note 创建原子类型对的三角矩阵，预计算所有可能的相互作用
     */
    precalculate(const ScoringFunction& sf, fl v=max_fl, fl factor=32){
        m_factor = factor;
        m_cutoff_sqr = sqr(sf.get_cutoff());
        m_max_cutoff_sqr = sqr(sf.get_max_cutoff());
        m_n = sz(m_factor * m_max_cutoff_sqr) + 3; // 确保索引范围安全
        
        // 创建原子类型对的三角矩阵
        triangular_matrix<precalculate_element> data(num_atom_types(sf.get_atom_typing()), precalculate_element(m_n, m_factor));

        VINA_CHECK(m_factor > epsilon_fl);
        VINA_CHECK(sz(m_max_cutoff_sqr * m_factor) + 1 < m_n);
        VINA_CHECK(m_max_cutoff_sqr * m_factor + 1 < m_n);

        flv rs = calculate_rs(); // 计算距离数组

        // 为每对原子类型计算相互作用能量
        VINA_FOR(t1, data.dim())
        {
            VINA_RANGE(t2, t1, data.dim())
            {
                precalculate_element &p = data(t1, t2);
                // 初始化smooth[].first (能量值)
                VINA_FOR_IN(i, p.smooth)
                {
                    p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));
                }

                // 初始化其余部分（导数和快速查找表）
                p.init_from_smooth_fst(rs);
            }
        }

        m_data = data;
    };
    
    /**
     * @brief 快速能量评估
     * @param type_pair_index 原子类型对索引
     * @param r2 距离的平方
     * @return 相互作用能量
     */
    fl eval_fast(sz type_pair_index, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(type_pair_index).eval_fast(r2);
    };
    
    /**
     * @brief 评估能量及其导数
     * @param type_pair_index 原子类型对索引
     * @param r2 距离的平方
     * @return pair<能量值, 导数>
     */
    pr eval_deriv(sz type_pair_index, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(type_pair_index).eval_deriv(r2);
    };
    
    /**
     * @brief 获取原子类型对的索引（容错版本）
     * @param t1 第一个原子类型
     * @param t2 第二个原子类型
     * @return 类型对索引
     */
    sz index_permissive(sz t1, sz t2) const { return m_data.index_permissive(t1, t2); }
    
    /**
     * @brief 获取截断距离的平方
     * @return 截断距离的平方
     */
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    
    /**
     * @brief 获取最大截断距离的平方
     * @return 最大截断距离的平方
     */
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    
    /**
     * @brief 扩展所有原子类型对的能量函数
     * @param left 向左扩展的距离
     * @param right 向右扩展的距离
     */
    void widen(fl left, fl right){
        flv rs = calculate_rs();
        VINA_FOR(t1, m_data.dim())
        VINA_RANGE(t2, t1, m_data.dim())
        m_data(t1, t2).widen(rs, left, right);
    };
};

/**
 * @class precalculate_byatom
 * @brief 基于具体原子的预计算类
 * 
 * 为分子模型中的每对具体原子预计算相互作用能量，用于更精确的局部优化
 */
struct precalculate_byatom
{
private:
    /**
     * @brief 计算距离数组
     * @return 距离值数组
     */
    flv calculate_rs() const{
        flv tmp(m_n, 0);
        VINA_FOR(i, m_n)
            tmp[i] = std::sqrt(i / m_factor);
        return tmp;
    };

    fl m_cutoff_sqr; ///< 截断距离的平方
    fl m_max_cutoff_sqr; ///< 最大截断距离的平方
    sz m_n; ///< 查找表大小
    fl m_factor; ///< 距离平方到索引的转换因子

    triangular_matrix<precalculate_element> m_data; ///< 存储所有原子对预计算数据的三角矩阵
public:
    /**
     * @brief 默认构造函数
     */
    precalculate_byatom() { }
    
    /**
     * @brief 构造函数，基于评分函数和分子模型初始化
     * @param sf 评分函数对象（应为连续函数以保证导数计算准确）
     * @param model 分子模型
     * @param v 能量上限值，默认为最大浮点数
     * @param factor 距离平方到索引的转换因子，默认为32
     * @note 为模型中每对原子预计算相互作用，比基于类型的预计算更精确但内存消耗更大
     */
    precalculate_byatom(const ScoringFunction &sf, const model &model, fl v=max_fl, fl factor=32){
        m_factor = factor;
        m_cutoff_sqr = sqr(sf.get_cutoff());
        m_max_cutoff_sqr = sqr(sf.get_max_cutoff());
        m_n = sz(m_factor * m_max_cutoff_sqr) + 3;
        
        sz n_atoms = model.num_atoms();
        atomv atoms = model.get_atoms();
        // 创建原子对的三角矩阵
        triangular_matrix<precalculate_element> data(n_atoms, precalculate_element(m_n, m_factor));

        VINA_CHECK(m_factor > epsilon_fl);
        VINA_CHECK(sz(m_max_cutoff_sqr * m_factor) + 1 < m_n);
        VINA_CHECK(m_max_cutoff_sqr * m_factor + 1 < m_n);

        flv rs = calculate_rs();

        // 为每对原子计算相互作用能量
        VINA_FOR(i, data.dim())
        {
            VINA_RANGE(j, i, data.dim())
            {
                precalculate_element &p = data(i, j);
                // 初始化smooth[].first (能量值)
                VINA_FOR_IN(k, p.smooth)
                {
                    p.smooth[k].first = (std::min)(v, sf.eval(atoms[i], atoms[j], rs[k]));
                }

                // 初始化其余部分
                p.init_from_smooth_fst(rs);
            }
        }
        m_data = data;
    };
    
    /**
     * @brief 快速能量评估
     * @param i 第一个原子索引
     * @param j 第二个原子索引
     * @param r2 距离的平方
     * @return 相互作用能量
     */
    fl eval_fast(sz i, sz j, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(i, j).eval_fast(r2);
    };
    
    /**
     * @brief 评估能量及其导数
     * @param i 第一个原子索引
     * @param j 第二个原子索引
     * @param r2 距离的平方
     * @return pair<能量值, 导数>
     */
    pr eval_deriv(sz i, sz j, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(i, j).eval_deriv(r2);
    };
    
    /**
     * @brief 获取原子对的索引（容错版本）
     * @param t1 第一个原子索引
     * @param t2 第二个原子索引
     * @return 原子对索引
     */
    sz index_permissive(sz t1, sz t2) const { return m_data.index_permissive(t1, t2); }
    
    /**
     * @brief 获取截断距离的平方
     * @return 截断距离的平方
     */
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    
    /**
     * @brief 获取最大截断距离的平方
     * @return 最大截断距离的平方
     */
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    
    /**
     * @brief 获取转换因子
     * @return 距离平方到索引的转换因子
     */
    sz get_factor() const { return m_factor; }
    
    /**
     * @brief 扩展所有原子对的能量函数
     * @param left 向左扩展的距离
     * @param right 向右扩展的距离
     */
    void widen(fl left, fl right){
        flv rs = calculate_rs();
        VINA_FOR(t1, m_data.dim())
        VINA_RANGE(t2, t1, m_data.dim())
        m_data(t1, t2).widen(rs, left, right);
    };
};

#endif
