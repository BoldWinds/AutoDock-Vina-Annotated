/**
 * @file conf.h
 * @brief AutoDock Vina分子构象配置和操作函数定义
*/

/*

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

#ifndef VINA_CONF_H
#define VINA_CONF_H

#include <boost/ptr_container/ptr_vector.hpp> // typedef output_container

#include "quaternion.h"
#include "random.h"

/**
 * @brief 缩放因子结构体，用于控制分子运动的步长
 */
struct scale {
	fl position;		///< 位置变化的缩放因子
	fl orientation;		///< 取向变化的缩放因子
	fl torsion;			 ///< 扭转角变化的缩放因子
	scale(fl position_, fl orientation_, fl torsion_) : position(position_), orientation(orientation_), torsion(torsion_) {}
};

/**
 * @brief 构象尺寸描述结构体，记录配体和柔性残基的扭转角数量
 * @note 用于初始化构象和变化结构体
 */
struct conf_size {
	szv ligands;	///< 每个配体的扭转角数量向量
	szv flex;		///< 每个柔性残基（受体）的扭转角数量向量

	/**
     * @brief 计算总的自由度数量
     * @return 总自由度数（配体扭转角 + 柔性残基扭转角 + 配体刚体运动6*N）
     */
	sz num_degrees_of_freedom() const {
		return sum(ligands) + sum(flex) + 6 * ligands.size();
	}
};

/**
 * @brief 将所有扭转角设置为0
 * @param torsions 扭转角向量
 */
inline void torsions_set_to_null(flv& torsions) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = 0;
}


/**
 * @brief 按比例增量更新扭转角并进行角度归一化
 * @param torsions 目标扭转角向量
 * @param c 增量向量
 * @param factor 缩放因子
 * @note 新的扭转角会被归一化到[-π, π]范围内
 */
inline void torsions_increment(flv& torsions, const flv& c, fl factor) { // new torsions are normalized
	VINA_FOR_IN(i, torsions) {
		torsions[i] += normalized_angle(factor * c[i]);
		normalize_angle(torsions[i]);
	}
}

/**
 * @brief 随机化扭转角
 * @param torsions 扭转角向量
 * @param generator 随机数生成器
 * @note 每个扭转角被设为[-π, π]范围内的随机值
 */
inline void torsions_randomize(flv& torsions, rng& generator) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = random_fl(-pi, pi, generator);
}

/**
 * @brief 检查两组扭转角是否过于接近
 * @param torsions1 第一组扭转角
 * @param torsions2 第二组扭转角
 * @param cutoff 角度阈值
 * @return 如果所有对应角度差都小于阈值则返回true
 */
inline bool torsions_too_close(const flv& torsions1, const flv& torsions2, fl cutoff) {
	assert(torsions1.size() == torsions2.size());
	VINA_FOR_IN(i, torsions1)
		if(std::abs(normalized_angle(torsions1[i] - torsions2[i])) > cutoff) 
			return false;
	return true;
}

/**
 * @brief 生成新的扭转角构象
 * @param torsions 目标扭转角向量
 * @param spread 随机扰动范围
 * @param rp 参考构象使用概率
 * @param rs 参考构象指针（可为空）
 * @param generator 随机数生成器
 * @note 对每个扭转角，以rp概率使用参考值，否则在当前值基础上加随机扰动
 */
inline void torsions_generate(flv& torsions, fl spread, fl rp, const flv* rs, rng& generator) {
	assert(!rs || rs->size() == torsions.size()); // if present, rs should be the same size as torsions
	VINA_FOR_IN(i, torsions)
		if(rs && random_fl(0, 1, generator) < rp)
			torsions[i] = (*rs)[i];
		else
			torsions[i] += random_fl(-spread, spread, generator);
}

/**
 * @brief 刚体变化结构体，描述分子刚体部分的位置和取向变化
 */
struct rigid_change {
	vec position;
	vec orientation;
	rigid_change() : position(0, 0, 0), orientation(0, 0, 0) {}
	void print() const {
		::print(position);
		::print(orientation);
	}
};


/**
 * @brief 刚体构象结构体，描述分子刚体部分的位置和取向
 */
struct rigid_conf {
	vec position;
	qt orientation;
	rigid_conf() : position(0, 0, 0), orientation(qt_identity) {}
	/**
     * @brief 重置为初始状态
     */
	void set_to_null() {
		position = zero_vec;
		orientation = qt_identity;
	}
	/**
     * @brief 按比例增量更新刚体构象
     * @param c 变化量
     * @param factor 缩放因子
     */
	void increment(const rigid_change& c, fl factor) {
		position += factor * c.position;
		vec rotation; rotation = factor * c.orientation;
		quaternion_increment(orientation, rotation); // orientation does not get normalized; tests show rounding errors growing very slowly
	}
	/**
     * @brief 在指定盒子内随机化位置和取向
     * @param corner1 盒子一角
     * @param corner2 盒子对角
     * @param generator 随机数生成器
     */
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		position = random_in_box(corner1, corner2, generator);
		orientation = random_orientation(generator);
	}
	/**
     * @brief 检查与另一构象是否过于接近
     * @param c 比较的构象
     * @param position_cutoff 位置阈值
     * @param orientation_cutoff 取向阈值
     * @return 如果位置和取向差异都小于阈值则返回true
     */
	bool too_close(const rigid_conf& c, fl position_cutoff, fl orientation_cutoff) const {
		if(vec_distance_sqr(position, c.position) > sqr(position_cutoff)) return false;
		if(sqr(quaternion_difference(orientation, c.orientation)) > sqr(orientation_cutoff)) return false;
		return true;
	}
	/**
     * @brief 随机扰动位置
     * @param spread 扰动范围
     * @param generator 随机数生成器
     */
	void mutate_position(fl spread, rng& generator) {
		position += spread * random_inside_sphere(generator);
	}
	/**
     * @brief 随机扰动取向
     * @param spread 扰动范围
     * @param generator 随机数生成器
     */
	void mutate_orientation(fl spread, rng& generator) {
		vec tmp; tmp = spread * random_inside_sphere(generator);
		quaternion_increment(orientation, tmp);
	}
	/**
     * @brief 生成新的刚体构象
     * @param position_spread 位置扰动范围
     * @param orientation_spread 取向扰动范围
     * @param rp 参考构象使用概率
     * @param rs 参考构象指针
     * @param generator 随机数生成器
	 * @note 在BFGS算法中，往往会根据一定概率选择最佳构象(rs)或者进行局部扰动
     */
	void generate(fl position_spread, fl orientation_spread, fl rp, const rigid_conf* rs, rng& generator) {
		if(rs && random_fl(0, 1, generator) < rp)
			position = rs->position;
		else
			mutate_position(position_spread, generator);
		if(rs && random_fl(0, 1, generator) < rp)
			orientation = rs->orientation;
		else
			mutate_orientation(orientation_spread, generator);
	}
	/**
     * @brief 将刚体变换应用到坐标集合
     * @param in 输入坐标，常量，往往是以分子质心为原点的标准构象
     * @param out 输出坐标，是输入坐标经过变换后的实际空间坐标
     * @param begin 起始索引
     * @param end 结束索引
     * @note 执行旋转+平移变换：out = R*in + t；似乎没有代码调用这个函数？
     */
	void apply(const vecv& in, vecv& out, sz begin, sz end) const {
		assert(in.size() == out.size());
		const mat m = quaternion_to_r3(orientation);
		VINA_RANGE(i, begin, end)
			out[i] = m * in[i] + position;
	}
	void print() const {
		::print(position);
		::print(orientation);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & position;
		ar & orientation;
	}
};

/**
 * @brief 配体变化结构体，包含刚体变化和扭转角变化
 */
struct ligand_change {
	rigid_change rigid;
	flv torsions;
	void print() const {
		rigid.print();
		printnl(torsions);
	}
};

/**
 * @brief 配体构象结构体，完整描述配体的构象状态
 */
struct ligand_conf {
	rigid_conf rigid;
	flv torsions;
	/**
     * @brief 重置为初始状态
     */
	void set_to_null() {
		rigid.set_to_null();
		torsions_set_to_null(torsions);
	}
	/**
     * @brief 按比例增量更新配体构象
     * @param c 变化量
     * @param factor 缩放因子
     */
	void increment(const ligand_change& c, fl factor) {
		rigid.increment(c.rigid, factor);
		torsions_increment(torsions, c.torsions, factor);
	}
	/**
     * @brief 随机化配体构象
     * @param corner1 位置盒子一角
     * @param corner2 位置盒子对角
     * @param generator 随机数生成器
     */
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		rigid.randomize(corner1, corner2, generator);
		torsions_randomize(torsions, generator);
	}
	/**
     * @brief 打印构象信息
     */
	void print() const {
		rigid.print();
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & rigid;
		ar & torsions;
	}
};

/**
 * @brief 残基变化结构体，仅包含扭转角变化
 */
struct residue_change {
	flv torsions;
	void print() const {
		printnl(torsions);
	}
};

/**
 * @brief 残基构象结构体，描述柔性残基的构象状态
 */
struct residue_conf {
	flv torsions;
	/**
     * @brief 重置为初始状态
     */
	void set_to_null() {
		torsions_set_to_null(torsions);
	}
	/**
     * @brief 按比例增量更新残基构象
     * @param c 变化量
     * @param factor 缩放因子
     */
	void increment(const residue_change& c, fl factor) {
		torsions_increment(torsions, c.torsions, factor);
	}
	/**
     * @brief 随机化残基构象
     * @param generator 随机数生成器
     */
	void randomize(rng& generator) {
		torsions_randomize(torsions, generator);
	}
	/**
     * @brief 打印构象信息
     */
	void print() const {
		printnl(torsions);
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & torsions;
	}
};

/**
 * @brief 整体变化结构体，包含所有配体和柔性残基的变化
 */
struct change {
	std::vector<ligand_change> ligands;
	std::vector<residue_change> flex;
	/**
     * @brief 构造函数，根据尺寸规格初始化
     * @param s 构象尺寸描述
     */
	change(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0);
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0);
	}
	/**
     * @brief 按索引访问变化值（只读）
     * @param index 线性索引
     * @return 对应的变化值
     * @note 索引顺序：配体1(位置3+取向3+扭转角N1), 配体2(...), 柔性残基1(扭转角M1), ...
     */
	fl operator()(sz index) const { // returns by value
		VINA_FOR_IN(i, ligands) {
			const ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			const residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return 0; // shouldn't happen, placating the compiler
	}
	/**
     * @brief 按索引访问变化值（可写）
     * @param index 线性索引
     * @return 对应变化值的引用
     */
	fl& operator()(sz index) {
		VINA_FOR_IN(i, ligands) {
			ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return ligands[0].rigid.position[0]; // shouldn't happen, placating the compiler
	}
	/**
     * @brief 获取总浮点数个数
     * @return 所有变化参数的总数
     */
	sz num_floats() const {
		sz tmp = 0;
		VINA_FOR_IN(i, ligands)
			tmp += 6 + ligands[i].torsions.size();
		VINA_FOR_IN(i, flex)
			tmp += flex[i].torsions.size();
		return tmp;
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
};

/**
 * @brief 完整构象结构体，描述系统的完整构象状态
 */
struct conf {
	std::vector<ligand_conf> ligands;	///< 配体构象向量
	std::vector<residue_conf> flex;		///< 柔性残基构象向量

	conf() {}
	    /**
     * @brief 构造函数，根据尺寸规格初始化
     * @param s 构象尺寸描述
     */
	conf(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0); // FIXME?
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0); // FIXME?
	}
	/**
     * @brief 重置为初始状态
     */
	void set_to_null() {
		VINA_FOR_IN(i, ligands)
			ligands[i].set_to_null();
		VINA_FOR_IN(i, flex)
			flex[i].set_to_null();
	}
	/**
     * @brief 按比例增量更新构象
     * @param c 变化量
     * @param factor 缩放因子
     * @note 扭转角会被归一化，但取向不会
     */
	void increment(const change& c, fl factor) { // torsions get normalized, orientations do not
		VINA_FOR_IN(i, ligands)
			ligands[i].increment(c.ligands[i], factor);
		VINA_FOR_IN(i, flex)
			flex[i]   .increment(c.flex[i],    factor);
	}
	/**
     * @brief 检查内部扭转角是否过于接近
     * @param c 比较的构象
     * @param torsions_cutoff 扭转角阈值
     * @return 如果所有扭转角差异都小于阈值则返回true
     */
	bool internal_too_close(const conf& c, fl torsions_cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff))
				return false;
		return true;
	}
	/**
     * @brief 检查外部自由度是否过于接近
     * @param c 比较的构象
     * @param cutoff 各类阈值
     * @return 如果位置、取向和柔性残基扭转角差异都小于阈值则返回true
     */
	bool external_too_close(const conf& c, const scale& cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!ligands[i].rigid.too_close(c.ligands[i].rigid, cutoff.position, cutoff.orientation))
				return false;
		assert(flex.size() == c.flex.size());
		VINA_FOR_IN(i, flex)
			if(!torsions_too_close(flex[i].torsions, c.flex[i].torsions, cutoff.torsion))
				return false;
		return true;
	}
	/**
     * @brief 检查两构象是否过于接近
     * @param c 比较的构象
     * @param cutoff 各类阈值
     * @return 如果内部和外部自由度都过于接近则返回true
     */
	bool too_close(const conf& c, const scale& cutoff) const {
		return internal_too_close(c, cutoff.torsion) &&
			   external_too_close(c, cutoff); // a more efficient implementation is possible, probably
	}
	/**
     * @brief 生成内部构象（仅扭转角）
     * @param torsion_spread 扭转角扰动范围
     * @param rp 参考构象使用概率
     * @param rs 参考构象指针
     * @param generator 随机数生成器
     * @note 扭转角在此函数后不会被归一化
     */
	void generate_internal(fl torsion_spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position.assign(0);
			ligands[i].rigid.orientation = qt_identity;
			const flv* torsions_rs = rs ? (&rs->ligands[i].torsions) : NULL;
			torsions_generate(ligands[i].torsions, torsion_spread, rp, torsions_rs, generator);
		}
	}
	/**
     * @brief 生成外部构象（位置、取向、柔性残基扭转角）
     * @param spread 各类扰动范围
     * @param rp 参考构象使用概率
     * @param rs 参考构象指针
     * @param generator 随机数生成器
     * @note 扭转角在此函数后不会被归一化
     */
	void generate_external(const scale& spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			const rigid_conf* rigid_conf_rs = rs ? (&rs->ligands[i].rigid) : NULL;
			ligands[i].rigid.generate(spread.position, spread.orientation, rp, rigid_conf_rs, generator);
		}
		VINA_FOR_IN(i, flex) {
			const flv* torsions_rs = rs ? (&rs->flex[i].torsions) : NULL;
			torsions_generate(flex[i].torsions, spread.torsion, rp, torsions_rs, generator);
		}
	}
	/**
     * @brief 完全随机化构象
     * @param corner1 位置盒子一角
     * @param corner2 位置盒子对角
     * @param generator 随机数生成器
     */
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		VINA_FOR_IN(i, ligands)
			ligands[i].randomize(corner1, corner2, generator);
		VINA_FOR_IN(i, flex)
			flex[i].randomize(generator);
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & ligands;
		ar & flex;
	}
};

/**
 * @brief 输出结果结构体，包含构象和能量信息
 */
struct output_type {
	conf c;                ///< 构象
	fl e;                  ///< 结合能
	fl lb;                 ///< 下界能量
	fl ub;                 ///< 上界能量
	fl intra;              ///< 分子内能量
	fl inter;              ///< 分子间能量
	fl conf_independent;   ///< 构象无关能量
	fl unbound;            ///< 未结合状态能量
	fl total;              ///< 总能量
	vecv coords;           ///< 原子坐标

	output_type(const conf& c_, fl e_) : c(c_), e(e_) {}
	//output_type(const conf& c_, fl e_, fl intra_, fl conf_independent_) : c(c_), e(e_), intra(intra_), conf_independent(conf_independent_) {}
};

typedef boost::ptr_vector<output_type> output_container;	///< 输出结果容器类型

/**
 * @brief 比较运算符，用于结果排序
 * @param a 第一个结果
 * @param b 第二个结果
 * @return 根据结合能排序
 */
inline bool operator<(const output_type& a, const output_type& b) { // for sorting output_container
	return a.e < b.e;
}

#endif
