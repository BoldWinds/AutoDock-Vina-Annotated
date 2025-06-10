/**
 * @file conf_independent.h
 * @brief 构象无关输入参数和评分函数的定义
 * 
 * 该文件定义了AutoDock Vina中与分子构象无关的特征参数和相应的评分函数类。
 * 这些参数包括扭转角数量、重原子数量、疏水原子数量等分子固有属性。

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

#ifndef VINA_CONF_INDEPENDENT_H
#define VINA_CONF_INDEPENDENT_H

#include <stdlib.h>
#include "atom.h"
#include "common.h"

// 前向声明
struct model;

/**
 * @class conf_independent_inputs
 * @brief 构象无关输入参数类
 * 
 * 存储分子的构象无关特征，这些特征在分子对接过程中保持不变，
 * 用于计算构象无关的评分项。
 */
class conf_independent_inputs {
public:
    fl torsdof;                    ///< 来自PDBQT文件TORSDOF关键字的扭转自由度数量
    fl num_tors;                   ///< 可旋转键数量
    fl num_rotors;                 ///< 转子数量
    fl num_heavy_atoms;            ///< 重原子（非氢原子）数量
    fl num_hydrophobic_atoms;      ///< 疏水原子数量
    fl ligand_max_num_h_bonds;     ///< 配体最大氢键数量
    fl num_ligands;                ///< 配体数量
    fl ligand_lengths_sum;         ///< 配体长度总和

    /**
     * @brief 转换为浮点数向量
     * @return 包含所有参数的浮点数向量
     */
    operator flv() const;
    
    /**
     * @brief 默认构造函数
     * 初始化所有参数为默认值
     */
    conf_independent_inputs();
    
    /**
     * @brief 从模型构造
     * @param m 分子模型对象
     * 根据分子模型计算并初始化所有构象无关参数
     */
    conf_independent_inputs(const model& m);
    
    /**
     * @brief 获取参数名称列表
     * @return 包含所有参数名称的字符串向量
     */
    std::vector<std::string> get_names() const;
    
private:
    /**
     * @brief 计算原子的重原子邻居数量
     * @param m 分子模型
     * @param i 原子索引
     * @return 与该原子键合的重原子数量
     */
    unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const;
    
    /**
     * @brief 计算原子的转子数量
     * @param m 分子模型
     * @param i 原子索引
     * @return 该原子连接的可旋转键数量（连接到重配体原子）
     */
    unsigned atom_rotors(const model& m, const atom_index& i) const;
};

/**
 * @class ConfIndependent
 * @brief 构象无关评分函数基类
 * 
 * 定义了构象无关评分项的统一接口，所有具体的构象无关评分函数都继承自此基类。
 */
class ConfIndependent {
public:
    /**
     * @brief 虚析构函数
     * 确保派生类对象能够正确销毁
     */
    virtual ~ConfIndependent() { }
    
    /**
     * @brief 评估构象无关评分
     * @param in 构象无关输入参数
     * @param x 当前能量
     * @param i 迭代器（用于多参数情况）
     * @return 修正后的能量
     * @note 基类默认返回0，派生类需要重写此方法
     */
    virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) { return 0; };
};

/**
 * @class num_tors_sqr
 * @brief 扭转角数量平方项评分类
 * 
 * Vina评分函数中的一项，计算扭转角数量的平方对评分的贡献。
 * 用于惩罚具有过多可旋转键的分子。
 */
class num_tors_sqr : public ConfIndependent {
public:
    /**
     * @brief 默认构造函数
     */
    num_tors_sqr() { }
    
    /**
     * @brief 获取参数数量
     * @return 参数数量（固定为1）
     */
    sz size() const { return 1; }
    
    /**
     * @brief 计算扭转角平方项评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_tors_sqrt
 * @brief 扭转角数量开方项评分类
 * 
 * 计算扭转角数量的平方根对评分的贡献。
 */
class num_tors_sqrt : public ConfIndependent {
public:
    num_tors_sqrt() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算扭转角开方项评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_tors_div
 * @brief 扭转角数量除法项评分类
 * 
 * 计算基于扭转角数量的除法运算对评分的贡献。
 */
class num_tors_div : public ConfIndependent {
public:
    num_tors_div() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算扭转角除法项评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class ligand_length
 * @brief 配体长度评分类
 * 
 * 计算配体分子长度对评分的贡献。
 */
class ligand_length : public ConfIndependent {
public:
    ligand_length() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算配体长度评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_ligands
 * @brief 配体数量评分类
 * 
 * 计算配体数量对评分的贡献。
 */
class num_ligands : public ConfIndependent {
public:
    num_ligands() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算配体数量评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_heavy_atoms_div
 * @brief 重原子数量除法项评分类
 * 
 * 计算基于重原子数量的除法运算对评分的贡献。
 */
class num_heavy_atoms_div : public ConfIndependent {
public:
    num_heavy_atoms_div() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算重原子除法项评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_heavy_atoms
 * @brief 重原子数量评分类
 * 
 * 计算重原子数量对评分的贡献。
 */
class num_heavy_atoms : public ConfIndependent {
public:
    num_heavy_atoms() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算重原子数量评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class num_hydrophobic_atoms
 * @brief 疏水原子数量评分类
 * 
 * 计算疏水原子数量对评分的贡献。
 */
class num_hydrophobic_atoms : public ConfIndependent {
public:
    num_hydrophobic_atoms() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算疏水原子数量评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

/**
 * @class ad4_tors_add
 * @brief AutoDock4扭转角加法项评分类
 * 
 * 计算AutoDock4兼容的扭转角加法项对评分的贡献。
 */
class ad4_tors_add : public ConfIndependent {
public:
    ad4_tors_add() { }
    sz size() const { return 1; }
    
    /**
     * @brief 计算AD4扭转角加法项评分
     * @param in 构象无关输入参数
     * @param x 权重系数
     * @param i 迭代器
     * @return 评分值
     */
    fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i);
};

#endif
