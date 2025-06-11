/**
 * @file atom.h
 * @brief AutoDock Vina原子数据结构定义

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

#ifndef VINA_ATOM_H
#define VINA_ATOM_H

#include "atom_base.h"

/**
 * @struct atom_index
 * @brief 原子索引结构体
 * 
 * 用于标识和定位原子的索引信息，包含原子编号和网格位置标志
 */
struct atom_index {
	sz i;           ///< 原子索引编号
	bool in_grid;   ///< 是否位于计算网格内的标
	/**
     * @brief 默认构造函数
     * @note 初始化索引为最大值，网格标志为false
     */
	atom_index() : i(max_sz), in_grid(false) {}
	/**
     * @brief 参数化构造函数
     * @param i_ 原子索引编号
     * @param in_grid_ 是否在网格内的标志
     */
	atom_index(sz i_, bool in_grid_) : i(i_), in_grid(in_grid_) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & i;
		ar & in_grid;
	}
};

/**
 * @brief 原子索引相等性比较运算符
 * @param i 第一个原子索引
 * @param j 第二个原子索引
 * @return 如果两个索引的编号和网格标志都相同则返回true
 */
inline bool operator==(const atom_index& i, const atom_index& j) {
	return i.i == j.i && i.in_grid == j.in_grid;
}

/**
 * @struct bond
 * @brief 化学键结构体
 * 
 * 表示两个原子之间的化学键信息，包含连接原子、键长和旋转性
 */
struct bond {
	atom_index connected_atom_index;  ///< 连接的原子索引
	fl length;                        ///< 键长度
	bool rotatable;                   ///< 是否为可旋转键
	/**
     * @brief 默认构造函数
     * @note 初始化键长为0，可旋转性为false
     */
	bond() : length(0), rotatable(false) {}

	/**
     * @brief 参数化构造函数
     * @param connected_atom_index_ 连接的原子索引
     * @param length_ 键长度
     * @param rotatable_ 是否可旋转
     */
	bond(const atom_index& connected_atom_index_, fl length_, bool rotatable_) : connected_atom_index(connected_atom_index_), length(length_), rotatable(rotatable_) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & connected_atom_index;
		ar & length;
		ar & rotatable;
	}
};

/**
 * @struct atom
 * @brief 原子结构体
 * 
 * 继承自atom_base的完整原子表示，包含空间坐标和化学键信息
 * @note 用于分子对接计算中的原子建模
 */
struct atom : public atom_base {
	vec coords;                    ///< 原子的三维空间坐标
	std::vector<bond> bonds;       ///< 该原子参与的所有化学键列表
	/**
     * @brief 默认构造函数
     * @note 坐标初始化为最大向量值
     */
    atom() : coords(max_vec) {}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<atom_base>(*this);
		ar & coords;
		ar & bonds;
	}
};

typedef std::vector<atom> atomv;	///< 原子向量类型定义，用于存储多个原子

#endif
