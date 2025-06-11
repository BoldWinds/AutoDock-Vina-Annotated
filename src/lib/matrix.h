/**
 * @file matrix.h
 * @brief 矩阵模板类实现文件
 * @details 提供通用矩阵、三角矩阵和严格三角矩阵的模板类实现，用于AutoDock Vina的数值计算
 * 
 * Copyright (c) 2006-2010, The Scripps Research Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * Author: Dr. Oleg Trott <ot14@columbia.edu>, 
 *         The Olson Lab, 
 *         The Scripps Research Institute
 */

#ifndef VINA_MATRIX_H
#define VINA_MATRIX_H

#include <vector>
#include "triangular_matrix_index.h"

/**
 * @brief 矩阵访问操作符宏定义
 * @details 为减少重复代码，定义统一的矩阵元素访问操作符
 * @note 提供一维和二维索引访问方式，包括const和非const版本
 */
#define VINA_MATRIX_DEFINE_OPERATORS \
    const T& operator()(sz i) const { return m_data[i]; } \
          T& operator()(sz i)       { return m_data[i]; } \
    const T& operator()(sz i, sz j) const { return m_data[index(i, j)]; } \
          T& operator()(sz i, sz j)       { return m_data[index(i, j)]; }

/**
 * @brief 通用矩阵模板类
 * @tparam T 矩阵元素类型
 * @details 实现二维矩阵的存储和操作，采用列主序存储方式
 */
template<typename T>
class matrix {
    std::vector<T> m_data;  ///< 矩阵数据存储向量
    sz m_i, m_j;           ///< 矩阵维度：m_i为行数，m_j为列数
public:
    /**
     * @brief 计算二维索引对应的一维数组索引
     * @param i 行索引
     * @param j 列索引
     * @return sz 一维数组中对应的索引位置
     * @note 采用列主序存储：index = i + m_i*j
     */
    sz index(sz i, sz j) const {
        assert(j < m_j);  // 列索引有效性检查
        assert(i < m_i);  // 行索引有效性检查
        return i + m_i*j; // 列主序索引计算
    }
    
    /**
     * @brief 默认构造函数
     * @details 创建空矩阵（0x0）
     */
    matrix() : m_i(0), m_j(0) {}
    
    /**
     * @brief 带参数构造函数
     * @param i 矩阵行数
     * @param j 矩阵列数
     * @param filler_val 初始填充值
     */
    matrix(sz i, sz j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {}
    
    /**
     * @brief 重新调整矩阵大小
     * @param m 新的行数
     * @param n 新的列数
     * @param filler_val 新增元素的填充值
     * @note 新尺寸必须大于等于原尺寸，保留原有数据
     */
    void resize(sz m, sz n, const T& filler_val) {
        if(m == dim_1() && n == dim_2()) return; // 尺寸未变，直接返回
        VINA_CHECK(m >= dim_1());  // 新行数必须不小于原行数
        VINA_CHECK(n >= dim_2());  // 新列数必须不小于原列数
        std::vector<T> tmp(m*n, filler_val);  // 创建新的存储空间
        // 复制原有数据到新空间
        VINA_FOR(i, m_i)
            VINA_FOR(j, m_j)
                tmp[i+m*j] = (*this)(i, j);
        m_data = tmp;  // 更新数据存储
        m_i = m;       // 更新行数
        m_j = n;       // 更新列数
    }
    
    /**
     * @brief 追加另一个矩阵到当前矩阵
     * @param x 要追加的矩阵
     * @param filler_val 填充值
     * @details 将矩阵x追加到当前矩阵的右下角，形成更大的矩阵
     */
    void append(const matrix<T>& x, const T& filler_val) {
        sz m = dim_1();  // 当前矩阵行数
        sz n = dim_2();  // 当前矩阵列数
        resize(m + x.dim_1(), n + x.dim_2(), filler_val);  // 扩展矩阵尺寸
        // 复制追加矩阵的数据到右下角
        VINA_FOR(i, x.dim_1())
            VINA_FOR(j, x.dim_2())
                (*this)(i+m, j+n) = x(i, j);
    }
    
    VINA_MATRIX_DEFINE_OPERATORS  ///< 展开矩阵访问操作符宏
    
    /**
     * @brief 获取矩阵行数
     * @return sz 矩阵行数
     */
    sz dim_1() const { return m_i; }
    
    /**
     * @brief 获取矩阵列数
     * @return sz 矩阵列数
     */
    sz dim_2() const { return m_j; }
};

/**
 * @brief 三角矩阵模板类
 * @tparam T 矩阵元素类型
 * @details 实现对称三角矩阵的存储和操作，只存储上三角部分（包括对角线）
 */
template<typename T>
class triangular_matrix {
    std::vector<T> m_data;  ///< 三角矩阵数据存储向量
    sz m_dim;              ///< 矩阵维度（n×n矩阵）
public:
    /**
     * @brief 计算三角矩阵的索引
     * @param i 行索引
     * @param j 列索引
     * @return sz 一维数组中对应的索引位置
     */
    sz index(sz i, sz j) const { return triangular_matrix_index(m_dim, i, j); }
    
    /**
     * @brief 宽松索引计算（自动处理i>j的情况）
     * @param i 行索引
     * @param j 列索引
     * @return sz 一维数组中对应的索引位置
     * @note 如果i>j，自动交换i和j的位置
     */
    sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }
    
    /**
     * @brief 默认构造函数
     */
    triangular_matrix() : m_dim(0) {}
    
    /**
     * @brief 带参数构造函数
     * @param n 矩阵维度
     * @param filler_val 初始填充值
     * @note 存储空间大小为n*(n+1)/2
     */
    triangular_matrix(sz n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {} 
    
    VINA_MATRIX_DEFINE_OPERATORS  ///< 展开矩阵访问操作符宏
    
    /**
     * @brief 获取矩阵维度
     * @return sz 矩阵维度
     */
    sz dim() const { return m_dim; }
};

/**
 * @brief 严格三角矩阵模板类
 * @tparam T 矩阵元素类型
 * @details 实现严格上三角矩阵的存储和操作，不包括对角线元素
 */
template<typename T>
class strictly_triangular_matrix {
    std::vector<T> m_data;  ///< 严格三角矩阵数据存储向量
    sz m_dim;              ///< 矩阵维度（n×n矩阵）
public:
    /**
     * @brief 计算严格三角矩阵的索引
     * @param i 行索引
     * @param j 列索引
     * @return sz 一维数组中对应的索引位置
     * @note 要求i < j，即只存储上三角部分（不含对角线）
     */
    sz index(sz i, sz j) const {
        assert(j < m_dim);  // 列索引有效性检查
        assert(i < j);      // 必须是上三角位置
        assert(j >= 1);     // 列索引至少为1
        return i + j*(j-1)/2;  // 严格三角矩阵索引计算
    }
    
    /**
     * @brief 宽松索引计算（自动处理i>j的情况）
     * @param i 行索引
     * @param j 列索引
     * @return sz 一维数组中对应的索引位置
     */
    sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }
    
    /**
     * @brief 默认构造函数
     */
    strictly_triangular_matrix() : m_dim(0) {}
    
    /**
     * @brief 带参数构造函数
     * @param n 矩阵维度
     * @param filler_val 初始填充值
     * @note 存储空间大小为n*(n-1)/2
     */
    strictly_triangular_matrix(sz n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {}
    
    /**
     * @brief 重新调整矩阵大小
     * @param n 新的矩阵维度
     * @param filler_val 新增元素的填充值
     * @note 新维度必须大于原维度，保留原有数据
     */
    void resize(sz n, const T& filler_val) {
        if(n == m_dim) return;  // 尺寸未变，直接返回
        VINA_CHECK(n > m_dim);  // 新维度必须大于原维度
        m_dim = n;
        m_data.resize(n*(n-1)/2, filler_val);  // 扩展存储空间并保留原数据
    }
    
    /**
     * @brief 追加另一个严格三角矩阵
     * @param m 要追加的严格三角矩阵
     * @param filler_val 填充值
	 * 
     * @example 示例说明：
     * 假设当前严格三角矩阵（3×3）:
     *     0   1   2
     * 0   -   1   2
     * 1   -   -   3
     * 2   -   -   -
     * 
     * 要追加的严格三角矩阵（2×2）:
     *     0   1
     * 0   -   4
     * 1   -   -
     * 
     * 执行过程：
     * 1. 调用resize(3+2, 0)将当前矩阵扩展到5×5，新位置用填充值0填充
     * 2. 将追加矩阵的数据复制到位置(i+3, j+3)
     * 
     * 最终结果（5×5）:
     *       0   1   2   3   4
     *   0   -   1   2   0   0
     *   1   -   -   3   0   0
     *   2   -   -   -   0   0
     *   3   -   -   -   -   4
     *   4   -   -   -   -   -
     * 
     * 内部存储变化:
     * 原始: [1, 2, 3]
     * 扩展后: [1, 2, 3, 0, 0, 0, 0, 0, 0, 4]
     */
    void append(const strictly_triangular_matrix<T>& m, const T& filler_val) { 
        sz n = dim();  // 当前矩阵维度
        resize(n + m.dim(), filler_val);  // 扩展矩阵维度
        // 复制追加矩阵的数据
        VINA_FOR(i, m.dim())
            VINA_RANGE(j, i+1, m.dim())
                (*this)(i+n, j+n) = m(i, j);
    }
    
    /**
     * @brief 追加矩形矩阵和三角矩阵的组合
     * @param rectangular 矩形矩阵部分
     * @param triangular 三角矩阵部分
     * @details 将矩形矩阵和三角矩阵组合追加到当前矩阵
	 * 
     * @example 示例说明：
     * 假设当前严格三角矩阵（3×3）:
     *     0   1   2
     * 0   -   1   2
     * 1   -   -   3
     * 2   -   -   -
     * 
     * 矩形矩阵（3×2）:
     *     0   1
     * 0   5   6
     * 1   7   8
     * 2   9  10
     * 
     * 严格三角矩阵（2×2）:
     *     0   1
     * 0   -  11
     * 1   -   -
     * 
     * 执行过程：
     * 1. 先调用append(triangular, 5)将当前矩阵扩展到5×5
     * 2. 将矩形矩阵数据填入位置(i, n+j)，其中n=3
     * 
     * 最终结果（5×5）:
     *       0   1   2   3   4
     *   0   -   1   2   5   6
     *   1   -   -   3   7   8
     *   2   -   -   -   9  10
     *   3   -   -   -   -  11
     *   4   -   -   -   -   -
     */
    void append(const matrix<T>& rectangular, const strictly_triangular_matrix<T>& triangular) {
        VINA_CHECK(dim() == rectangular.dim_1());      // 行数必须匹配
        VINA_CHECK(rectangular.dim_2() == triangular.dim());  // 列数必须匹配

        // 处理特殊情况：矩形矩阵为空的情况
        if(rectangular.dim_2() == 0) return;
        if(rectangular.dim_1() == 0) {
            (*this) = triangular;  // 直接赋值三角矩阵
            return;
        }
        const T& filler_val = rectangular(0, 0);  // 使用矩形矩阵元素作为填充值

        sz n = dim();  // 当前矩阵维度
        append(triangular, filler_val);  // 先追加三角矩阵
        // 然后添加矩形矩阵数据
        VINA_FOR(i, rectangular.dim_1())
            VINA_FOR(j, rectangular.dim_2())
                (*this)(i, n + j) = rectangular(i, j);
    }
    
    VINA_MATRIX_DEFINE_OPERATORS  ///< 展开矩阵访问操作符宏
    
    /**
     * @brief 获取矩阵维度
     * @return sz 矩阵维度
     */
    sz dim() const { return m_dim; }
};

#undef VINA_MATRIX_DEFINE_OPERATORS  ///< 取消宏定义

#endif
