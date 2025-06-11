/**
 * @file atom_constants.h
 * @brief AutoDock Vina原子类型和属性常量定义文件
 * 
 * 定义了三种评分函数(AutoDock4, X-Score, DrugScore-CSD)使用的原子类型常量，
 * 原子属性参数，以及相关的类型转换和判断函数。

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

#ifndef VINA_ATOM_CONSTANTS_H
#define VINA_ATOM_CONSTANTS_H

#include "common.h"

// 基于SY_TYPE_*的元素类型定义，包括氢原子
const sz EL_TYPE_H    =  0;
const sz EL_TYPE_C    =  1;
const sz EL_TYPE_N    =  2;
const sz EL_TYPE_O    =  3;
const sz EL_TYPE_S    =  4;
const sz EL_TYPE_P    =  5;
const sz EL_TYPE_F    =  6;
const sz EL_TYPE_Cl   =  7;
const sz EL_TYPE_Br   =  8;
const sz EL_TYPE_I    =  9;
const sz EL_TYPE_Si   = 10; // Silicon
const sz EL_TYPE_At   = 11; // Astatine
const sz EL_TYPE_Met  = 12;	///< 金属元素类型
const sz EL_TYPE_Dummy= 13;	///< 虚拟原子类型
const sz EL_TYPE_SIZE = 14;	///< 元素类型总数

// AutoDock4原子类型定义
const sz AD_TYPE_C    =  0;	///< 碳原子(非芳香性)
const sz AD_TYPE_A    =  1;	///< 芳香碳原子
const sz AD_TYPE_N    =  2;	///< 氮原子(非氢键)
const sz AD_TYPE_O    =  3;	///< 氧原子(非氢键)
const sz AD_TYPE_P    =  4;
const sz AD_TYPE_S    =  5;
const sz AD_TYPE_H    =  6;	///< 非极性氢原子
const sz AD_TYPE_F    =  7;
const sz AD_TYPE_I    =  8;
const sz AD_TYPE_NA   =  9;	///< 氮原子(氢键受体)
const sz AD_TYPE_OA   = 10;	///< 氧原子(氢键受体)
const sz AD_TYPE_SA   = 11;	///< 硫原子(氢键受体)
const sz AD_TYPE_HD   = 12;	///< 氢原子(氢键供体)
const sz AD_TYPE_Mg   = 13;
const sz AD_TYPE_Mn   = 14;
const sz AD_TYPE_Zn   = 15;
const sz AD_TYPE_Ca   = 16;
const sz AD_TYPE_Fe   = 17;
const sz AD_TYPE_Cl   = 18;
const sz AD_TYPE_Br   = 19;
const sz AD_TYPE_Si   = 20; // Silicon
const sz AD_TYPE_At   = 21; // Astatine
const sz AD_TYPE_G0   = 22;	///< 环闭合虚拟原子G0
const sz AD_TYPE_G1   = 23;	///< 环闭合虚拟原子G1
const sz AD_TYPE_G2   = 24;	///< 环闭合虚拟原子G2
const sz AD_TYPE_G3   = 25;	///< 环闭合虚拟原子G3
const sz AD_TYPE_CG0  = 26;	///< 环碳原子CG0
const sz AD_TYPE_CG1  = 27;	///< 环碳原子CG1
const sz AD_TYPE_CG2  = 28;	///< 环碳原子CG2
const sz AD_TYPE_CG3  = 29;	///< 环碳原子CG3
const sz AD_TYPE_W    = 30;	///< 水合配体虚拟原子
const sz AD_TYPE_SIZE = 31;	///< AutoDock原子类型总数

// X-Score评分函数原子类型定义
const sz XS_TYPE_C_H   =  0;  ///< 疏水碳原子
const sz XS_TYPE_C_P   =  1;  ///< 极性碳原子
const sz XS_TYPE_N_P   =  2;  ///< 极性氮原子
const sz XS_TYPE_N_D   =  3;  ///< 氢键供体氮原子
const sz XS_TYPE_N_A   =  4;  ///< 氢键受体氮原子
const sz XS_TYPE_N_DA  =  5;  ///< 氢键供体/受体氮原子
const sz XS_TYPE_O_P   =  6;  ///< 极性氧原子
const sz XS_TYPE_O_D   =  7;  ///< 氢键供体氧原子
const sz XS_TYPE_O_A   =  8;  ///< 氢键受体氧原子
const sz XS_TYPE_O_DA  =  9;  ///< 氢键供体/受体氧原子
const sz XS_TYPE_S_P   = 10;  ///< 极性硫原子
const sz XS_TYPE_P_P   = 11;  ///< 极性磷原子
const sz XS_TYPE_F_H   = 12;  ///< 疏水氟原子
const sz XS_TYPE_Cl_H  = 13;  ///< 疏水氯原子
const sz XS_TYPE_Br_H  = 14;  ///< 疏水溴原子
const sz XS_TYPE_I_H   = 15;  ///< 疏水碘原子
const sz XS_TYPE_Si    = 16;  ///< 硅原子
const sz XS_TYPE_At    = 17;  ///< 砹原子
const sz XS_TYPE_Met_D = 18;  ///< 金属氢键供体
const sz XS_TYPE_C_H_CG0 = 19; ///< 环闭合疏水碳原子CG0
const sz XS_TYPE_C_P_CG0 = 20; ///< 环闭合极性碳原子CG0
const sz XS_TYPE_G0      = 21; ///< 环闭合虚拟原子G0
const sz XS_TYPE_C_H_CG1 = 22; ///< 环闭合疏水碳原子CG1
const sz XS_TYPE_C_P_CG1 = 23; ///< 环闭合极性碳原子CG1
const sz XS_TYPE_G1      = 24; ///< 环闭合虚拟原子G1
const sz XS_TYPE_C_H_CG2 = 25; ///< 环闭合疏水碳原子CG2
const sz XS_TYPE_C_P_CG2 = 26; ///< 环闭合极性碳原子CG2
const sz XS_TYPE_G2      = 27; ///< 环闭合虚拟原子G2
const sz XS_TYPE_C_H_CG3 = 28; ///< 环闭合疏水碳原子CG3
const sz XS_TYPE_C_P_CG3 = 29; ///< 环闭合极性碳原子CG3
const sz XS_TYPE_G3      = 30; ///< 环闭合虚拟原子G3
const sz XS_TYPE_W       = 31; ///< 水合配体虚拟原子
const sz XS_TYPE_SIZE    = 32; ///< X-Score原子类型总数

// DrugScore-CSD评分函数原子类型定义
const sz SY_TYPE_C_3   =  0;  ///< sp3杂化碳原子
const sz SY_TYPE_C_2   =  1;  ///< sp2杂化碳原子
const sz SY_TYPE_C_ar  =  2;  ///< 芳香碳原子
const sz SY_TYPE_C_cat =  3;  ///< 正离子碳原子
const sz SY_TYPE_N_3   =  4;  ///< sp3杂化氮原子
const sz SY_TYPE_N_ar  =  5;  ///< 芳香氮原子
const sz SY_TYPE_N_am  =  6;  ///< 酰胺氮原子
const sz SY_TYPE_N_pl3 =  7;  ///< 平面三配位氮原子
const sz SY_TYPE_O_3   =  8;  ///< sp3杂化氧原子
const sz SY_TYPE_O_2   =  9;  ///< sp2杂化氧原子
const sz SY_TYPE_O_co2 = 10;  ///< 羧基氧原子
const sz SY_TYPE_S     = 11;  ///< 硫原子
const sz SY_TYPE_P     = 12;  ///< 磷原子
const sz SY_TYPE_F     = 13;  ///< 氟原子
const sz SY_TYPE_Cl    = 14;  ///< 氯原子
const sz SY_TYPE_Br    = 15;  ///< 溴原子
const sz SY_TYPE_I     = 16;  ///< 碘原子
const sz SY_TYPE_Met   = 17;  ///< 金属原子
const sz SY_TYPE_SIZE  = 18;  ///< DrugScore原子类型总数

/**
 * @brief 原子种类属性结构体
 * 
 * 包含原子的物理化学属性参数，用于能量计算
 */
struct atom_kind {
    std::string name;        ///< 原子名称
    fl radius;              ///< 范德华半径(Å)
    fl depth;               ///< 范德华势阱深度(kcal/mol)
    fl hb_depth;            ///< 氢键势阱深度，负值表示受体，正值表示供体
    fl hb_radius;           ///< 氢键相互作用半径(Å)
    fl solvation;           ///< 溶剂化参数
    fl volume;              ///< 原子体积(Å³)
    fl covalent_radius;     ///< 共价半径(Å)，来源于维基百科
};

/**
 * @brief 原子种类数据数组
 * 
 * 基于编辑过的AD4_parameters.data文件生成，共价半径来自维基百科
 * @note 数组索引对应AD_TYPE_*常量值
 */
const atom_kind atom_kind_data[] = {
    // name, radius, depth, hb_depth, hb_r, solvation, volume, covalent radius
	{ "C", 2.00000, 0.15000, 0.0, 0.0,   -0.00143,   33.51030,   0.77}, //  0
	{ "A", 2.00000, 0.15000, 0.0, 0.0,   -0.00052,   33.51030,   0.77}, //  1
	{ "N", 1.75000, 0.16000, 0.0, 0.0,   -0.00162,   22.44930,   0.75}, //  2
	{ "O", 1.60000, 0.20000, 0.0, 0.0,   -0.00251,   17.15730,   0.73}, //  3
	{ "P", 2.10000, 0.20000, 0.0, 0.0,   -0.00110,   38.79240,   1.06}, //  4
	{ "S", 2.00000, 0.20000, 0.0, 0.0,   -0.00214,   33.51030,   1.02}, //  5
	{ "H", 1.00000, 0.02000, 0.0, 0.0,    0.00051,    0.00000,   0.37}, //  6
	{ "F", 1.54500, 0.08000, 0.0, 0.0,   -0.00110,   15.44800,   0.71}, //  7
	{ "I", 2.36000, 0.55000, 0.0, 0.0,   -0.00110,   55.05850,   1.33}, //  8
	{"NA", 1.75000, 0.16000,-5.0, 1.9,   -0.00162,   22.44930,   0.75}, //  9
	{"OA", 1.60000, 0.20000,-5.0, 1.9,   -0.00251,   17.15730,   0.73}, // 10
	{"SA", 2.00000, 0.20000,-1.0, 2.5,   -0.00214,   33.51030,   1.02}, // 11
	{"HD", 1.00000, 0.02000, 1.0, 0.0,    0.00051,    0.00000,   0.37}, // 12
	{"Mg", 0.65000, 0.87500, 0.0, 0.0,   -0.00110,    1.56000,   1.30}, // 13
	{"Mn", 0.65000, 0.87500, 0.0, 0.0,   -0.00110,    2.14000,   1.39}, // 14
	{"Zn", 0.74000, 0.55000, 0.0, 0.0,   -0.00110,    1.70000,   1.31}, // 15
	{"Ca", 0.99000, 0.55000, 0.0, 0.0,   -0.00110,    2.77000,   1.74}, // 16
	{"Fe", 0.65000, 0.01000, 0.0, 0.0,   -0.00110,    1.84000,   1.25}, // 17
	{"Cl", 2.04500, 0.27600, 0.0, 0.0,   -0.00110,   35.82350,   0.99}, // 18
	{"Br", 2.16500, 0.38900, 0.0, 0.0,   -0.00110,   42.56610,   1.14}, // 19
	{"Si", 2.30000, 0.20000, 0.0, 0.0,   -0.00143,   50.96500,   1.11}, // 20
	{"At", 2.40000, 0.55000, 0.0, 0.0,   -0.00110,   57.90580,   1.44}, // 21
	{"G0", 0.00000, 0.00000, 0.0, 0.0,    0.00000,    0.00000,   0.77}, // 22
	{"G1", 0.00000, 0.00000, 0.0, 0.0,    0.00000,    0.00000,   0.77}, // 23
	{"G2", 0.00000, 0.00000, 0.0, 0.0,    0.00000,    0.00000,   0.77}, // 24
	{"G3", 0.00000, 0.00000, 0.0, 0.0,    0.00000,    0.00000,   0.77}, // 25
	{"CG0",2.00000, 0.15000, 0.0, 0.0,   -0.00143,   33.51030,   0.77}, // 26
	{"CG1",2.00000, 0.15000, 0.0, 0.0,   -0.00143,   33.51030,   0.77}, // 27
	{"CG2",2.00000, 0.15000, 0.0, 0.0,   -0.00143,   33.51030,   0.77}, // 28
	{"CG3",2.00000, 0.15000, 0.0, 0.0,   -0.00143,   33.51030,   0.77}, // 29
	{"W",  0.00000, 0.00000, 0.0, 0.0,    0.00000,    0.00000,   0.00}  // 30
};

const fl metal_solvation_parameter = -0.00110;  ///< 金属原子通用溶剂化参数

const fl metal_covalent_radius = 1.75; // for metals not on the list // FIXME this info should be moved to non_ad_metals

const sz atom_kinds_size =  sizeof(atom_kind_data) / sizeof(const atom_kind);	///< 原子种类数据数组大小

/**
 * @brief 原子等价性结构体
 * 
 * 定义某些原子类型的等价替换关系
 */
struct atom_equivalence {
    std::string name;  ///< 原始原子名称
    std::string to;    ///< 等价替换的原子名称
};

/**
 * @brief 原子等价性数据数组
 */
const atom_equivalence atom_equivalence_data[] = {
    {"Se",  "S"}  ///< 硒原子等价为硫原子
};

const sz atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(const atom_equivalence);

/**
 * @brief 氢键受体种类结构体
 */
struct acceptor_kind {
	sz ad_type;	///< AutoDock原子类型
	fl radius;	///< 最适氢键长度(Å)
	fl depth;	///< 氢键势阱深度(kcal/mol)
};

/**
 * @brief 氢键受体种类数据数组
 */
const acceptor_kind acceptor_kind_data[] = { // ad_type, optimal length, depth
	{AD_TYPE_NA, 1.9, 5.0},
	{AD_TYPE_OA, 1.9, 5.0},
	{AD_TYPE_SA, 2.5, 1.0}
};

const sz acceptor_kinds_size = sizeof(acceptor_kind_data) / sizeof(acceptor_kind);

/**
 * @brief 判断是否为氢原子
 * @param ad AutoDock原子类型
 * @return 如果是氢原子返回true，否则返回false
 */
inline bool ad_is_hydrogen(sz ad) {
	return ad == AD_TYPE_H || ad == AD_TYPE_HD;
}

/**
 * @brief 判断是否为杂原子(非C/H原子)
 * @param ad AutoDock原子类型
 * @return 如果是杂原子返回true，否则返回false
 * @note 对于ad >= AD_TYPE_SIZE的情况返回false
 */
inline bool ad_is_heteroatom(sz ad) { // returns false for ad >= AD_TYPE_SIZE
	return ad != AD_TYPE_A && ad != AD_TYPE_C  && 
		   ad != AD_TYPE_H && ad != AD_TYPE_HD && 
		   ad < AD_TYPE_SIZE;
}

/**
 * @brief AutoDock原子类型转换为元素类型
 * @param t AutoDock原子类型
 * @return 对应的元素类型
 */
inline sz ad_type_to_el_type(sz t) {
	switch(t) {
		case AD_TYPE_C    : return EL_TYPE_C;
		case AD_TYPE_A    : return EL_TYPE_C;
		case AD_TYPE_N    : return EL_TYPE_N;
		case AD_TYPE_O    : return EL_TYPE_O;
		case AD_TYPE_P    : return EL_TYPE_P;
		case AD_TYPE_S    : return EL_TYPE_S;
		case AD_TYPE_H    : return EL_TYPE_H;
		case AD_TYPE_F    : return EL_TYPE_F;
		case AD_TYPE_I    : return EL_TYPE_I;
		case AD_TYPE_NA   : return EL_TYPE_N;
		case AD_TYPE_OA   : return EL_TYPE_O;
		case AD_TYPE_SA   : return EL_TYPE_S;
		case AD_TYPE_HD   : return EL_TYPE_H;
		case AD_TYPE_Mg   : return EL_TYPE_Met;
		case AD_TYPE_Mn   : return EL_TYPE_Met;
		case AD_TYPE_Zn   : return EL_TYPE_Met;
		case AD_TYPE_Ca   : return EL_TYPE_Met;
		case AD_TYPE_Fe   : return EL_TYPE_Met;
		case AD_TYPE_Cl   : return EL_TYPE_Cl;
		case AD_TYPE_Br   : return EL_TYPE_Br;
		case AD_TYPE_Si   : return EL_TYPE_Si;
		case AD_TYPE_At   : return EL_TYPE_At;
		case AD_TYPE_CG0  : return EL_TYPE_C;
		case AD_TYPE_CG1  : return EL_TYPE_C;
		case AD_TYPE_CG2  : return EL_TYPE_C;
		case AD_TYPE_CG3  : return EL_TYPE_C;
		case AD_TYPE_G0   : return EL_TYPE_Dummy;
		case AD_TYPE_G1   : return EL_TYPE_Dummy;
		case AD_TYPE_G2   : return EL_TYPE_Dummy;
		case AD_TYPE_G3   : return EL_TYPE_Dummy;
		case AD_TYPE_W    : return EL_TYPE_Dummy;
		case AD_TYPE_SIZE : return EL_TYPE_SIZE;
		default: VINA_CHECK(false);
	}
	return EL_TYPE_SIZE; // to placate the compiler in case of warnings - it should never get here though
}

/**
 * @brief X-Score评分函数范德华半径数组
 * @note 数组索引对应XS_TYPE_*常量值
 */
const fl xs_vdw_radii[] = {
	1.9, // C_H
	1.9, // C_P
	1.8, // N_P
	1.8, // N_D
	1.8, // N_A
	1.8, // N_DA
	1.7, // O_P
	1.7, // O_D
	1.7, // O_A
	1.7, // O_DA
	2.0, // S_P
	2.1, // P_P
	1.5, // F_H
	1.8, // Cl_H
	2.0, // Br_H
	2.2, // I_H
    2.2, // Si
    2.3, // At
	1.2, // Met_D
    1.9, // C_H_CG0
    1.9, // C_P_CG0
    1.9, // C_H_CG1
    1.9, // C_P_CG1
    1.9, // C_H_CG2
    1.9, // C_P_CG2
    1.9, // C_H_CG3
    1.9, // C_P_CG3
    0.0, // G0
    0.0, // G1
    0.0, // G2
    0.0, // G3
    0.0  // W
};

/**
 * @brief Vinardo评分函数范德华半径数组
 * @note 数组索引对应XS_TYPE_*常量值，Vinardo是Vina的改进版本
 */
const fl xs_vinardo_vdw_radii[] = {
	2.0, // C_H
	2.0, // C_P
	1.7, // N_P
	1.7, // N_D
	1.7, // N_A
	1.7, // N_DA
	1.6, // O_P
	1.6, // O_D
	1.6, // O_A
	1.6, // O_DA
	2.0, // S_P
	2.1, // P_P
	1.5, // F_H
	1.8, // Cl_H
	2.0, // Br_H
	2.2, // I_H
	2.2, // Si
	2.3, // At
	1.2, // Met_D
	2.0, // C_H_CG0
	2.0, // C_P_CG0
	2.0, // C_H_CG1
	2.0, // C_P_CG1
	2.0, // C_H_CG2
	2.0, // C_P_CG2
	2.0, // C_H_CG3
	2.0, // C_P_CG3
	0.0, // G0
	0.0, // G1
	0.0, // G2
	0.0, // G3
	0.0	 // W
};

/**
 * @brief 获取X-Score评分函数原子半径
 * @param t X-Score原子类型
 * @return 原子的范德华半径(Å)
 */
inline fl xs_radius(sz t) {
	const sz n = sizeof(xs_vdw_radii) / sizeof(const fl);
	assert(n == XS_TYPE_SIZE);
	assert(t < n);
	return xs_vdw_radii[t];
}

/**
 * @brief 获取Vinardo评分函数原子半径
 * @param t X-Score原子类型
 * @return 原子的范德华半径(Å)
 */
inline fl xs_vinardo_radius(sz t) {
	const sz n = sizeof(xs_vdw_radii) / sizeof(const fl);
	assert(n == XS_TYPE_SIZE);
	assert(t < n);
	return xs_vinardo_vdw_radii[t];
}

/**
 * @brief 非AutoDock金属原子名称数组
 * @note 包含AutoDock4不直接支持但需要特殊处理的金属原子
 */
const std::string non_ad_metal_names[] = { // expand as necessary
	"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

/**
 * @brief 判断是否为非AutoDock金属原子
 * @param name 原子名称
 * @return 如果是非AutoDock金属原子返回true，否则返回false
 */
inline bool is_non_ad_metal_name(const std::string& name) {
	const sz s = sizeof(non_ad_metal_names) / sizeof(const std::string);
	VINA_FOR(i, s)
		if(non_ad_metal_names[i] == name)
			return true;
	return false;
}

/**
 * @brief 判断X-Score原子类型是否为疏水性
 * @param xs X-Score原子类型
 * @return 如果是疏水原子返回true，否则返回false
 */
inline bool xs_is_hydrophobic(sz xs) {
	return xs == XS_TYPE_C_H || 
		   xs == XS_TYPE_F_H ||
		   xs == XS_TYPE_Cl_H ||
		   xs == XS_TYPE_Br_H || 
		   xs == XS_TYPE_I_H;
}

/**
 * @brief 判断X-Score原子类型是否为氢键受体
 * @param xs X-Score原子类型
 * @return 如果是氢键受体返回true，否则返回false
 */
inline bool xs_is_acceptor(sz xs) {
	return xs == XS_TYPE_N_A ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_A ||
		   xs == XS_TYPE_O_DA;
}

/**
 * @brief 判断X-Score原子类型是否为氢键供体
 * @param xs X-Score原子类型
 * @return 如果是氢键供体返回true，否则返回false
 */
inline bool xs_is_donor(sz xs) {
	return xs == XS_TYPE_N_D ||
		   xs == XS_TYPE_N_DA ||
		   xs == XS_TYPE_O_D ||
		   xs == XS_TYPE_O_DA ||
		   xs == XS_TYPE_Met_D;
}

/**
 * @brief 判断两个原子是否能形成供体-受体氢键
 * @param t1 第一个原子的X-Score类型
 * @param t2 第二个原子的X-Score类型
 * @return 如果t1是供体且t2是受体返回true，否则返回false
 */
inline bool xs_donor_acceptor(sz t1, sz t2) {
	return xs_is_donor(t1) && xs_is_acceptor(t2);
}

/**
 * @brief 判断两个原子是否可能形成氢键
 * @param t1 第一个原子的X-Score类型
 * @param t2 第二个原子的X-Score类型
 * @return 如果两个原子可能形成氢键返回true，否则返回false
 * @note 检查双向的供体-受体关系
 */
inline bool xs_h_bond_possible(sz t1, sz t2) {
	return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

/**
 * @brief 获取AutoDock原子类型的属性信息
 * @param i AutoDock原子类型索引
 * @return 对应的原子属性结构体引用
 */
inline const atom_kind& ad_type_property(sz i) {
	assert(AD_TYPE_SIZE == atom_kinds_size);
    assert(i < atom_kinds_size);
    return atom_kind_data[i];
}

/**
 * @brief 字符串原子名称转换为AutoDock原子类型
 * @param name 原子名称字符串
 * @return 对应的AutoDock原子类型，找不到时返回AD_TYPE_SIZE
 * @note 不抛出异常，因为未知金属原子并非异常情况
 */
inline sz string_to_ad_type(const std::string& name) { // returns AD_TYPE_SIZE if not found (no exceptions thrown, because metals unknown to AD4 are not exceptional)
    VINA_FOR(i, atom_kinds_size)
		if(atom_kind_data[i].name == name)
			return i;
	VINA_FOR(i, atom_equivalences_size)
		if(atom_equivalence_data[i].name == name)
			return string_to_ad_type(atom_equivalence_data[i].to);
    return AD_TYPE_SIZE;
}

/**
 * @brief 获取所有原子类型中的最大共价半径
 * @return 最大共价半径值(Å)
 */
inline fl max_covalent_radius() {
	fl tmp = 0;
	VINA_FOR(i, atom_kinds_size)
		if(atom_kind_data[i].covalent_radius > tmp)
			tmp = atom_kind_data[i].covalent_radius;
	return tmp;
}

#endif
