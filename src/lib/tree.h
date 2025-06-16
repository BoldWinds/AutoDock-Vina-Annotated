/**
 * @file tree.h
 * @brief AutoDock Vina分子树结构定义，用于表示和操作分子的层次结构

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

#ifndef VINA_TREE_H
#define VINA_TREE_H

#include "conf.h"
#include "atom.h"

/**
 * @brief 坐标系框架，提供原点和方向信息
 * @details 定义了分子片段的空间坐标系，包含原点位置和旋转方向
 */
struct frame {
	/**
     * @brief 构造函数，使用指定原点初始化框架
     * @param origin_ 框架原点位置
     */
	frame(const vec& origin_) : origin(origin_), orientation_q(qt_identity), orientation_m(quaternion_to_r3(qt_identity)) {}
	/**
     * @brief 将局部坐标转换为实验室坐标系
     * @param local_coords 局部坐标
     * @return 实验室坐标系中的位置
     */
	vec local_to_lab(const vec& local_coords) const {
		vec tmp;
		tmp = origin + orientation_m*local_coords; 
		return tmp;
	}
	/**
     * @brief 将局部方向向量转换为实验室坐标系
     * @param local_direction 局部方向向量
     * @return 实验室坐标系中的方向向量
     */
	vec local_to_lab_direction(const vec& local_direction) const {
		vec tmp;
		tmp = orientation_m * local_direction;
		return tmp;
	}
	/**
     * @brief 获取框架方向四元数
	 * @return 实验室坐标系中的方向四元数
     */
	const qt& orientation() const { return orientation_q; }
	/**
     * @brief 获取frame原点位置
     * @return frame原点位置
     */
	const vec& get_origin() const { return origin; }
protected:
	vec origin; ///< 框架原点位置
	/**
     * @brief 设置框架方向
     * @param q 方向四元数（不进行归一化）
     */
	void set_orientation(const qt& q) { // does not normalize the orientation
		orientation_q = q;
		orientation_m = quaternion_to_r3(orientation_q);
	}
	mat orientation_m; ///< 方向旋转矩阵
	qt  orientation_q; ///< 方向四元数
};

/**
 * @brief 原子范围结构，定义原子序列的起始和结束位置
 */
struct atom_range {
    sz begin;	///< 起始原子索引
    sz end;		///< 结束原子索引
	atom_range(sz begin_, sz end_) : begin(begin_), end(end_) {}
	/**
     * @brief 使用变换函数修改原子范围
     * @tparam F 变换函数类型
     * @param f 变换函数
     */
	template<typename F>
	void transform(const F& f) {
		sz diff = end - begin;
		begin = f(begin);
		end   = begin + diff;
	}
};

/**
 * @brief 原子框架，结合了坐标系框架和原子范围
 * @details 管理一组原子的坐标变换和力计算
 */
struct atom_frame : public frame, public atom_range {
	/**
     * @brief 构造函数
     * @param origin_ 框架原点
     * @param begin_ 起始原子索引
     * @param end_ 结束原子索引
     */
	atom_frame(const vec& origin_, sz begin_, sz end_) : frame(origin_), atom_range(begin_, end_) {}
	/**
     * @brief 设置原子坐标
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     */
	void set_coords(const atomv& atoms, vecv& coords) const {
		VINA_RANGE(i, begin, end)
			coords[i] = local_to_lab(atoms[i].coords);
	}
	/**
     * @brief 计算力和力矩的总和
     * @param coords 原子坐标
     * @param forces 作用力
     * @return 力和力矩的对
     */
	vecp sum_force_and_torque(const vecv& coords, const vecv& forces) const {
		vecp tmp;
		tmp.first.assign(0);
		tmp.second.assign(0);
		VINA_RANGE(i, begin, end) {
			tmp.first  += forces[i]; 
			tmp.second += cross_product(coords[i] - origin, forces[i]);
		}
		return tmp;
	}
};

/**
 * @brief 刚体结构，表示不可变形的分子片段
 */
struct rigid_body : public atom_frame {
	/**
     * @brief 构造函数
     * @param origin_ 刚体原点
     * @param begin_ 起始原子索引
     * @param end_ 结束原子索引
     */
	rigid_body(const vec& origin_, sz begin_, sz end_) : atom_frame(origin_, begin_, end_) {}
	/**
     * @brief 设置刚体构象
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 刚体构象
     */
	void set_conf(const atomv& atoms, vecv& coords, const rigid_conf& c) {
		origin = c.position;
		set_orientation(c.orientation);
		set_coords(atoms, coords);
	}
	/**
     * @brief 计算扭转角数量（刚体无扭转角）
     * @param s 扭转角计数器
     */
	void count_torsions(sz& s) const {} // do nothing
	/**
     * @brief 设置导数（梯度信息）
     * @param force_torque 力和力矩
     * @param c 刚体变化量（输出）
     */
	void set_derivative(const vecp& force_torque, rigid_change& c) const {
		c.position     = force_torque.first;
		c.orientation  = force_torque.second;
	}
};

/**
 * @brief 轴框架，定义绕轴旋转的分子片段
 */
struct axis_frame : public atom_frame {
	/**
     * @brief 构造函数
     * @param origin_ 框架原点
     * @param begin_ 起始原子索引
     * @param end_ 结束原子索引
     * @param axis_root 轴根点
     */
	axis_frame(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : atom_frame(origin_, begin_, end_) {
		vec diff; diff = origin - axis_root;
		fl nrm = diff.norm();
		VINA_CHECK(nrm >= epsilon_fl);
		axis = (1/nrm) * diff;
	}
	/**
     * @brief 设置扭转角导数
     * @param force_torque 力和力矩
     * @param c 扭转角变化量（输出）
     */
	void set_derivative(const vecp& force_torque, fl& c) const {
		c = force_torque.second * axis; // 力矩在轴方向的投影
	}
protected:
	vec axis; ///< 旋转轴单位向量
};

/**
 * @brief 分子片段，支持相对于父框架的扭转
 */
struct segment : public axis_frame {
	/**
     * @brief 构造函数
     * @param origin_ 片段原点
     * @param begin_ 起始原子索引
     * @param end_ 结束原子索引
     * @param axis_root 轴根点
     * @param parent 父框架
     * @note 要求父框架方向为单位四元数
     */
	segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root, const frame& parent) : axis_frame(origin_, begin_, end_, axis_root) {
		VINA_CHECK(eq(parent.orientation(), qt_identity)); // the only initial parent orientation this c'tor supports
		relative_axis = axis;
		relative_origin = origin - parent.get_origin();
	}
	/**
     * @brief 设置片段构象
     * @param parent 父框架
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 扭转角迭代器
     */
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		const fl torsion = *c;
		++c;
		origin = parent.local_to_lab(relative_origin);
		axis = parent.local_to_lab_direction(relative_axis);
		qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();
		quaternion_normalize_approx(tmp); // normalization added in 1.1.2
		//quaternion_normalize(tmp); // normalization added in 1.1.2
		set_orientation(tmp);
		set_coords(atoms, coords);
	}
	/**
     * @brief 计算扭转角数量
     * @param s 扭转角计数器
     */
	void count_torsions(sz& s) const {
		++s;
	}
private:
	vec relative_axis;		///< 相对轴向量
	vec relative_origin;	///< 相对原点位置
};

/**
 * @brief 首个片段，作为分子树的根节点
 */
struct first_segment : public axis_frame {
	/**
     * @brief 从普通片段构造首个片段
     * @param s 源片段
     */
	first_segment(const segment& s) : axis_frame(s) {}
	/**
     * @brief 构造函数
     * @param origin_ 片段原点
     * @param begin_ 起始原子索引
     * @param end_ 结束原子索引
     * @param axis_root 轴根点
     */
	first_segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : axis_frame(origin_, begin_, end_, axis_root) {}
	/**
     * @brief 设置首个片段构象
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param torsion 扭转角
     */
	void set_conf(const atomv& atoms, vecv& coords, fl torsion) {
		set_orientation(angle_to_quaternion(axis, torsion));
		set_coords(atoms, coords);
	}
	/**
     * @brief 计算扭转角数量
     * @param s 扭转角计数器
     */
	void count_torsions(sz& s) const {
		++s;
	}
};

/**
 * @brief 批量设置分支构象
 * @tparam T 分支类型
 * @param b 分支向量
 * @param parent 父框架
 * @param atoms 原子数组
 * @param coords 坐标数组（输出）
 * @param c 构象参数迭代器
 */
template<typename T> // T == branch
void branches_set_conf(std::vector<T>& b, const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
	VINA_FOR_IN(i, b)
		b[i].set_conf(parent, atoms, coords, c);
}

/**
 * @brief 计算分支的导数并累加到输出
 * @tparam T 分支类型
 * @param b 分支向量
 * @param origin 参考原点
 * @param coords 原子坐标
 * @param forces 作用力
 * @param out 输出力和力矩（累加）
 * @param d 导数迭代器
 */
template<typename T> // T == branch
void branches_derivative(const std::vector<T>& b, const vec& origin, const vecv& coords, const vecv& forces, vecp& out, flv::iterator& d) { // adds to out
	VINA_FOR_IN(i, b) {
		vecp force_torque = b[i].derivative(coords, forces, d);
		out.first  += force_torque.first;
		vec r; r = b[i].node.get_origin() - origin;
		out.second += cross_product(r, force_torque.first) + force_torque.second;
	}
}

/**
 * @brief 分子树结构，递归定义分子的层次结构
 * @tparam T 节点类型（通常为segment）
 */
template<typename T> // T == segment
struct tree {
	T node;	 ///< 当前节点
	std::vector< tree<T> > children; ///< 子树向量
	/**
     * @brief 构造函数
     * @param node_ 节点对象
     */
	tree(const T& node_) : node(node_) {}
	/**
     * @brief 设置树的构象
     * @param parent 父框架
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 构象参数迭代器
     */
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		node.set_conf(parent, atoms, coords, c);
		branches_set_conf(children, node, atoms, coords, c);
	}
	/**
     * @brief 梯度计算
     * @param coords 原子坐标
     * @param forces 作用力
     * @param p 导数参数迭代器
     * @return 力和力矩对
     */
	vecp derivative(const vecv& coords, const vecv& forces, flv::iterator& p) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		fl& d = *p; // reference
		++p;
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		return force_torque;
	}
};

typedef tree<segment> branch;			///< 分支类型定义
typedef std::vector<branch> branches;	///< 分支向量类型定义

/**
 * @brief 异构树结构，支持不同类型的根节点
 * @tparam Node 根节点类型（first_segment或rigid_body）
 */
template<typename Node> // Node == first_segment || rigid_body
struct heterotree {
	Node node;			///< 根节点
	branches children;	///< 子分支
	/**
     * @brief 构造函数
     * @param node_ 根节点对象
     */
	heterotree(const Node& node_) : node(node_) {}
	/**
     * @brief 设置配体构象
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 配体构象
     */
	void set_conf(const atomv& atoms, vecv& coords, const ligand_conf& c) {
		node.set_conf(atoms, coords, c.rigid);
		flv::const_iterator p = c.torsions.begin();
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());
	}
	/**
     * @brief 设置残基构象
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 残基构象
     */
	void set_conf(const atomv& atoms, vecv& coords, const residue_conf& c) {
		flv::const_iterator p = c.torsions.begin();
		node.set_conf(atoms, coords, *p);
		++p;
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());
	}
	/**
     * @brief 梯度计算
     * @param coords 原子坐标
     * @param forces 作用力
     * @param c 配体变化量（输出）
     */
	void derivative(const vecv& coords, const vecv& forces, ligand_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, c.rigid);
		assert(p == c.torsions.end());
	}
	/**
     * @brief 计算残基梯度
     * @param coords 原子坐标
     * @param forces 作用力
     * @param c 残基变化量（输出）
     */
	void derivative(const vecv& coords, const vecv& forces, residue_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		fl& d = *p; // reference
		++p;
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		assert(p == c.torsions.end());
	}
};

/**
 * @brief 递归计算树结构的扭转角数量
 * @tparam T 树类型 (main_branch, branch, flexible_body)
 * @param t 树对象
 * @param s 扭转角计数器
 */
template<typename T> // T = main_branch, branch, flexible_body
void count_torsions(const T& t, sz& s) {
	t.node.count_torsions(s);
	VINA_FOR_IN(i, t.children)
		count_torsions(t.children[i], s);
}

typedef heterotree<rigid_body> flexible_body;	///< 柔性体类型定义
typedef heterotree<first_segment> main_branch;	///< 主分支类型定义

/**
 * @brief 可变向量，扩展std::vector以支持批量操作
 * @tparam T 元素类型（flexible_body或main_branch）
 */
template<typename T> // T == flexible_body || main_branch
struct vector_mutable : public std::vector<T> {
	/**
     * @brief 批量设置构象
     * @tparam C 构象类型
     * @param atoms 原子数组
     * @param coords 坐标数组（输出）
     * @param c 构象向量
     */
	template<typename C>
	void set_conf(const atomv& atoms, vecv& coords, const std::vector<C>& c) { // C == ligand_conf || residue_conf
		VINA_FOR_IN(i, (*this))
			(*this)[i].set_conf(atoms, coords, c[i]);
	}
	/**
     * @brief 计算所有元素的扭转角数量
     * @return 扭转角数量向量
     */
	szv count_torsions() const {
		szv tmp(this->size(), 0);
		VINA_FOR_IN(i, (*this))
			::count_torsions((*this)[i], tmp[i]);
		return tmp;
	}
	/**
     * @brief 梯度计算
     * @tparam C 变化量类型
     * @param coords 原子坐标
     * @param forces 作用力
     * @param c 变化量向量（输出）
     */
	template<typename C>
	void derivative(const vecv& coords, const vecv& forces, std::vector<C>& c) const { // C == ligand_change || residue_change
		VINA_FOR_IN(i, (*this))
			(*this)[i].derivative(coords, forces, c[i]);
	}
};

/**
 * @brief 递归变换树或异构树结构中的原子范围
 * @tparam T 树类型
 * @tparam F 变换函数类型
 * @param t 树对象
 * @param f 变换函数
 */
template<typename T, typename F> // tree or heterotree - like structure
void transform_ranges(T& t, const F& f) {
	t.node.transform(f);
	VINA_FOR_IN(i, t.children)
		transform_ranges(t.children[i], f);
}

#endif
