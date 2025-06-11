/**
 * @file atom_type.h
 * @brief 原子类型结构体定义文件
 * 
 * 定义了原子类型结构体，支持多种原子分类方案(EL、AD、XS、SY)，
 * 提供原子类型判断、属性查询和类型对索引计算等功能。

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

#ifndef VINA_ATOM_TYPE_H
#define VINA_ATOM_TYPE_H

#include "atom_constants.h"
#include "triangular_matrix_index.h"

/**
 * @brief 原子类型结构体
 * 
 * 存储原子在不同分类方案下的类型信息，支持四种原子分类系统：
 * - EL: 元素类型
 * - AD: AutoDock4原子类型  
 * - XS: X-Score原子类型
 * - SY: DrugScore-CSD原子类型
 */
struct atom_type {
	/**
     * @brief 原子分类方案枚举
     */
	enum t {EL, AD, XS, SY};
	sz el, ad, xs, sy;
	/**
     * @brief 默认构造函数
     * 
     * 将所有类型索引初始化为对应的SIZE值，表示未分配状态
     */
	atom_type() : el(EL_TYPE_SIZE), ad(AD_TYPE_SIZE), xs(XS_TYPE_SIZE), sy(SY_TYPE_SIZE) {}
	/**
     * @brief 获取指定分类方案下的原子类型索引
     * @param atom_typing_used 要使用的原子分类方案
     * @return 对应分类方案下的原子类型索引
     */
	sz get(t atom_typing_used) const {
		switch(atom_typing_used) {
			case EL: return el;
			case AD: return ad;
			case XS: return xs;
			case SY: return sy;
			default: assert(false); return max_sz;
		}
	}
	/**
     * @brief 判断是否为氢原子
     * @return 如果是氢原子返回true，否则返回false
     * @note 基于AutoDock原子类型进行判断
     */
	bool is_hydrogen() const {
		return ad_is_hydrogen(ad);
	}
	/**
     * @brief 判断是否为杂原子(非C/H原子)
     * @return 如果是杂原子返回true，否则返回false
     * @note 同时检查AutoDock类型和X-Score金属供体类型
     */
	bool is_heteroatom() const {
		return ad_is_heteroatom(ad) || xs == XS_TYPE_Met_D;
	}
	/**
     * @brief 判断原子类型是否可接受
     * @return 如果原子类型有效返回true，否则返回false
     * @note 有效条件：AutoDock类型已分配 OR 为X-Score金属供体
     */
	bool acceptable_type() const {
		return ad < AD_TYPE_SIZE || xs == XS_TYPE_Met_D;
	}
	/**
     * @brief 分配元素类型
     * 
     * 根据AutoDock类型自动分配对应的元素类型，
     * 特殊处理X-Score金属供体类型
     */
	void assign_el() {
		el = ad_type_to_el_type(ad);
		if(ad == AD_TYPE_SIZE && xs == XS_TYPE_Met_D)
			el = EL_TYPE_Met;
	}
	/**
     * @brief 判断两个原子是否为同一元素
     * @param a 另一个原子类型
     * @return 如果是同一元素返回true，否则返回false
     * @note 不区分不同金属或未分配类型
     */
	bool same_element(const atom_type& a) const { // does not distinguish metals or unassigned types
		return el == a.el;
	}
	/**
     * @brief 获取原子的共价半径
     * @return 原子的共价半径值(Å)
     * @note 优先使用AutoDock参数，特殊处理X-Score金属供体
     */
	fl covalent_radius() const {
		if(ad < AD_TYPE_SIZE)        return ad_type_property(ad).covalent_radius;
		else if(xs == XS_TYPE_Met_D) return metal_covalent_radius;
		VINA_CHECK(false);           return 0; // never happens - placating the compiler
	}
	/**
     * @brief 计算与另一个原子的最适共价键长度
     * @param x 另一个原子类型
     * @return 两个原子间的最适共价键长度(Å)
     * @note 简单地将两个原子的共价半径相加
     */
	fl optimal_covalent_bond_length(const atom_type& x) const {
		return covalent_radius() + x.covalent_radius();
	}
private:
	friend class boost::serialization::access;
	template<class Archive> 
	void serialize(Archive& ar, const unsigned version) {
		ar & el;
		ar & ad;
		ar & xs;
		ar & sy;
	}
};

/**
 * @brief 获取指定原子分类方案的类型总数
 * @param atom_typing_used 原子分类方案
 * @return 该分类方案下的原子类型总数
 */
inline sz num_atom_types(atom_type::t atom_typing_used) {
	switch(atom_typing_used) {
		case atom_type::EL: return EL_TYPE_SIZE;
		case atom_type::AD: return AD_TYPE_SIZE;
		case atom_type::XS: return XS_TYPE_SIZE;
		case atom_type::SY: return SY_TYPE_SIZE;
		default: assert(false); return max_sz;
	}
}

/**
 * @brief 计算原子对在三角矩阵中的索引
 * @param atom_typing_used 使用的原子分类方案
 * @param a 第一个原子类型
 * @param b 第二个原子类型
 * @return 原子对在三角矩阵中的索引
 * @throws 如果任一原子类型在指定方案下未分配则抛出错误
 * @note 获取索引，用于在预计算时查找原子间相互作用参数。
 */
inline sz get_type_pair_index(atom_type::t atom_typing_used, const atom_type& a, const atom_type& b) { // throws error if any arg is unassigned in the given typing scheme
	sz n = num_atom_types(atom_typing_used);

	sz i = a.get(atom_typing_used); VINA_CHECK(i < n);
	sz j = b.get(atom_typing_used); VINA_CHECK(j < n);

	if(i <= j) return triangular_matrix_index(n, i, j);
	else       return triangular_matrix_index(n, j, i);
}

#endif
