/**
 * @file conf.h
 * @brief AutoDock Vina分子构象配置和操作函数定义

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


/// @brief 缩放因子结构体，用于控制分子运动的步长
struct scale {
	fl position;		///< 位置变化的缩放因子
	fl orientation;		///< 取向变化的缩放因子
	fl torsion;			 ///< 扭转角变化的缩放因子
	scale(fl position_, fl orientation_, fl torsion_) : position(position_), orientation(orientation_), torsion(torsion_) {}
};

/// @brief 构象尺寸描述结构体，记录配体和柔性残基的扭转角数量，用于初始化构象和变化结构体
struct conf_size {
	szv ligands;	///< 每个配体的扭转角数量向量
	szv flex;		///< 每个柔性残基（受体）的扭转角数量向量

	/// @brief 计算总的自由度数量
	sz num_degrees_of_freedom() const {
		return sum(ligands) + sum(flex) + 6 * ligands.size();
	}
};

/// @brief 将所有扭转角设置为0
inline void torsions_set_to_null(flv& torsions) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = 0;
}


/// @brief 按比例增量更新扭转角并进行角度归一化
inline void torsions_increment(flv& torsions, const flv& c, fl factor) { // new torsions are normalized
	VINA_FOR_IN(i, torsions) {
		torsions[i] += normalized_angle(factor * c[i]);
		normalize_angle(torsions[i]);
	}
}

/// @brief 随机化扭转角，每个扭转角被设为[-π, π]范围内的随机值
inline void torsions_randomize(flv& torsions, rng& generator) {
	VINA_FOR_IN(i, torsions)
		torsions[i] = random_fl(-pi, pi, generator);
}

/// @brief 检查两组扭转角是否过于接近
inline bool torsions_too_close(const flv& torsions1, const flv& torsions2, fl cutoff) {
	assert(torsions1.size() == torsions2.size());
	VINA_FOR_IN(i, torsions1)
		if(std::abs(normalized_angle(torsions1[i] - torsions2[i])) > cutoff) 
			return false;
	return true;
}

/// @brief 生成新的扭转角构象：rp的概率直接使用参考构象rs，其他在当前值的基础上添加spread大小的扰动
/// @note 似乎目前没用
inline void torsions_generate(flv& torsions, fl spread, fl rp, const flv* rs, rng& generator) {
	assert(!rs || rs->size() == torsions.size()); // if present, rs should be the same size as torsions
	VINA_FOR_IN(i, torsions)
		if(rs && random_fl(0, 1, generator) < rp)
			torsions[i] = (*rs)[i];
		else
			torsions[i] += random_fl(-spread, spread, generator);
}

/// @brief 刚体变化结构体，描述分子刚体部分的位置和取向变化
struct rigid_change {
	vec position;
	vec orientation;
	rigid_change() : position(0, 0, 0), orientation(0, 0, 0) {}
	void print() const {
		::print(position);
		::print(orientation);
	}
};


/// @brief 刚体构象结构体，描述分子刚体部分的位置和取向
struct rigid_conf {
	vec position;
	qt orientation;
	rigid_conf() : position(0, 0, 0), orientation(qt_identity) {}
	/// @brief 重置为初始状态
	void set_to_null() {
		position = zero_vec;
		orientation = qt_identity;
	}
	/// @brief 按比例增量更新刚体构象
	void increment(const rigid_change& c, fl factor) {
		position += factor * c.position;
		vec rotation; rotation = factor * c.orientation;
		quaternion_increment(orientation, rotation); // orientation does not get normalized; tests show rounding errors growing very slowly
	}
	/// @brief 在指定盒子内随机化位置和取向
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		position = random_in_box(corner1, corner2, generator);
		orientation = random_orientation(generator);
	}
	/// @brief 检查与另一构象是否过于接近
	bool too_close(const rigid_conf& c, fl position_cutoff, fl orientation_cutoff) const {
		if(vec_distance_sqr(position, c.position) > sqr(position_cutoff)) return false;
		if(sqr(quaternion_difference(orientation, c.orientation)) > sqr(orientation_cutoff)) return false;
		return true;
	}
	/// @brief 随机扰动位置
	void mutate_position(fl spread, rng& generator) {
		position += spread * random_inside_sphere(generator);
	}
	/// @brief 随机扰动朝向
	void mutate_orientation(fl spread, rng& generator) {
		vec tmp; tmp = spread * random_inside_sphere(generator);
		quaternion_increment(orientation, tmp);
	}
	/// @brief 生成新的刚体构象
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
	/// @brief 将刚体变换应用到坐标集合
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

/// @brief 配体变化结构体，包含刚体变化和扭转角变化
struct ligand_change {
	rigid_change rigid;
	flv torsions;
	void print() const {
		rigid.print();
		printnl(torsions);
	}
};

/// @brief 配体构象结构体，完整描述配体的构象状态
struct ligand_conf {
	rigid_conf rigid;
	flv torsions;
	/// @brief 重置为初始状态
	void set_to_null() {
		rigid.set_to_null();
		torsions_set_to_null(torsions);
	}
	/// @brief 按比例增量更新配体构象
	void increment(const ligand_change& c, fl factor) {
		rigid.increment(c.rigid, factor);
		torsions_increment(torsions, c.torsions, factor);
	}
	/// @brief 随机化配体构象
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		rigid.randomize(corner1, corner2, generator);
		torsions_randomize(torsions, generator);
	}
	/// @brief 打印构象信息
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

/// @brief 残基变化结构体，仅包含扭转角变化
struct residue_change {
	flv torsions;
	void print() const {
		printnl(torsions);
	}
};

/// @brief 残基构象结构体，描述柔性残基的构象状态
struct residue_conf {
	flv torsions;
	/// @brief 重置为初始状态
	void set_to_null() {
		torsions_set_to_null(torsions);
	}
	/// @brief 按比例增量更新残基构象
	void increment(const residue_change& c, fl factor) {
		torsions_increment(torsions, c.torsions, factor);
	}
	/// @brief 随机化残基构象
	void randomize(rng& generator) {
		torsions_randomize(torsions, generator);
	}
	/// @brief 打印构象信息
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

/// @brief 整体变化结构体，包含所有配体和柔性残基的变化
struct change {
	std::vector<ligand_change> ligands;
	std::vector<residue_change> flex;
	/// @brief 构造函数，根据conf_size初始化
	change(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0);
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0);
	}
    /// @brief 按索引访问变化值
    /// @note 索引顺序：配体1(位置3+取向3+扭转角N1), 配体2(...), 柔性残基1(扭转角M1), ...
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
    /// @brief 按索引访问变化引用
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
	/// @brief 获取所有变化量的总数
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

/// @brief 完整构象结构体，描述系统的完整构象状态
struct conf {
	std::vector<ligand_conf> ligands;	///< 配体构象向量
	std::vector<residue_conf> flex;		///< 柔性残基构象向量

	conf() {}
	/// @brief 构造函数，根据尺寸规格初始化
	conf(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0); // FIXME?
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0); // FIXME?
	}
	/// @brief 重置为初始状态
	void set_to_null() {
		VINA_FOR_IN(i, ligands)
			ligands[i].set_to_null();
		VINA_FOR_IN(i, flex)
			flex[i].set_to_null();
	}
	/// @brief 按比例增量更新构象
	void increment(const change& c, fl factor) { // torsions get normalized, orientations do not
		VINA_FOR_IN(i, ligands)
			ligands[i].increment(c.ligands[i], factor);
		VINA_FOR_IN(i, flex)
			flex[i]   .increment(c.flex[i],    factor);
	}
	/// @brief 检查内部扭转角是否过于接近，仅在too_close中被调用
	bool internal_too_close(const conf& c, fl torsions_cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff))
				return false;
		return true;
	}
	/// @brief 检查外部自由度是否过于接近，仅在too_close中被调用
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
	/// 检查两构象是否过于接近
	bool too_close(const conf& c, const scale& cutoff) const {
		return internal_too_close(c, cutoff.torsion) &&
			   external_too_close(c, cutoff); // a more efficient implementation is possible, probably
	}
	/*	没有函数需要调用这部分
	/// @brief 生成内部构象（仅扭转角）
	void generate_internal(fl torsion_spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position.assign(0);
			ligands[i].rigid.orientation = qt_identity;
			const flv* torsions_rs = rs ? (&rs->ligands[i].torsions) : NULL;
			torsions_generate(ligands[i].torsions, torsion_spread, rp, torsions_rs, generator);
		}
	}
    /// @brief 生成外部构象（位置、取向、柔性残基扭转角）
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
	*/
	/// @brief 在盒子中随机化生成构象
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

/// @brief 输出结果结构体，包含构象和能量信息
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
};

typedef boost::ptr_vector<output_type> output_container;	///< 输出结果容器类型

/// @brief 比较两个输出结果的能量，用于排序
inline bool operator<(const output_type& a, const output_type& b) { // for sorting output_container
	return a.e < b.e;
}

#endif
