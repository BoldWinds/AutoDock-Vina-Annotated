/**
 * @file array3d.h
 * @brief 三维数组模板类的实现，提供安全的内存分配和访问功能

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

#ifndef VINA_ARRAY3D_H
#define VINA_ARRAY3D_H

#include <exception> // std::bad_alloc
#include "common.h"

/**
 * @brief 检查两个数值相乘是否会导致溢出
 * @param i 第一个乘数
 * @param j 第二个乘数
 * @return 乘积结果
 * @throws std::bad_alloc 当乘法溢出时抛出异常
 * @note 通过检查乘积是否小于任一乘数来检测溢出
 */
inline sz checked_multiply(sz i, sz j) {
	if(i == 0 || j == 0) return 0;
	const sz tmp = i * j;
	if(tmp < i || tmp < j || tmp / i != j)
		throw std::bad_alloc(); // can't alloc if the size makes sz wrap around
	return tmp;
}

/**
 * @brief 检查三个数值相乘是否会导致溢出
 * @param i 第一个乘数
 * @param j 第二个乘数
 * @param k 第三个乘数
 * @return 三个数的乘积
 * @throws std::bad_alloc 当乘法溢出时抛出异常
 * @note 递归调用二元checked_multiply函数实现三元检查
 */
inline sz checked_multiply(sz i, sz j, sz k) {
    return checked_multiply(checked_multiply(i, j), k);
}

/**
 * @brief 三维数组模板类
 * @tparam T 数组元素类型
 * @note 使用一维vector作为底层存储，通过索引计算实现三维访问
 */
template<typename T>
class array3d {
	sz m_i, m_j, m_k;           ///< 三个维度的大小
	std::vector<T> m_data;      ///< 底层一维数据存储

    // Boost序列化支持
    friend class boost::serialization::access;
    /**
     * @brief 序列化函数，用于数据持久化
     * @tparam Archive 归档类型
     * @param ar 归档对象
     * @param version 版本号（未使用）
     */
	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & m_i;
		ar & m_j;
		ar & m_k;
		ar & m_data;
	}
public:
    /**
     * @brief 默认构造函数，创建空的三维数组
     */
	array3d() : m_i(0), m_j(0), m_k(0) {}

    /**
     * @brief 构造指定大小的三维数组
     * @param i 第一维大小
     * @param j 第二维大小
     * @param k 第三维大小
     * @throws std::bad_alloc 当总大小溢出时抛出异常
     */
	array3d(sz i, sz j, sz k) : m_i(i), m_j(j), m_k(k), m_data(checked_multiply(i, j, k)) {}

    /**
     * @brief 获取第一维大小
     * @return 第一维的大小
     */
	sz dim0() const { return m_i; }

    /**
     * @brief 获取第二维大小
     * @return 第二维的大小
     */
	sz dim1() const { return m_j; }

    /**
     * @brief 获取第三维大小
     * @return 第三维的大小
     */
	sz dim2() const { return m_k; }

    /**
     * @brief 获取指定维度的大小
     * @param i 维度索引（0, 1, 2）
     * @return 对应维度的大小
     * @note 超出范围时会触发断言失败
     */
	sz dim(sz i) const {
		switch(i) {
			case 0: return m_i;
			case 1: return m_j;
			case 2: return m_k;
			default: assert(false); return 0; // to get rid of the warning
		}
	}

    /**
     * @brief 重新调整数组大小
     * @param i 新的第一维大小
     * @param j 新的第二维大小
     * @param k 新的第三维大小
     * @warning 原有数据将被清除
     */
	void resize(sz i, sz j, sz k) { // data is essentially garbled
		m_i = i;
		m_j = j;
		m_k = k;
		m_data.resize(checked_multiply(i, j, k));
	}

    /**
     * @brief 访问指定位置的元素（非const版本）
     * @param i 第一维索引
     * @param j 第二维索引
     * @param k 第三维索引
     * @return 元素的引用
     * @note 使用行优先存储顺序：index = i + m_i*(j + m_j*k)
     */
	T&       operator()(sz i, sz j, sz k)       { return m_data[i + m_i*(j + m_j*k)]; }

    /**
     * @brief 访问指定位置的元素（const版本）
     * @param i 第一维索引
     * @param j 第二维索引
     * @param k 第三维索引
     * @return 元素的常量引用
     * @note 使用行优先存储顺序：index = i + m_i*(j + m_j*k)
     */
	const T& operator()(sz i, sz j, sz k) const { return m_data[i + m_i*(j + m_j*k)]; }
};

#endif
