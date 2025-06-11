/**
 * @file conf_independent.cpp
 * @brief 构象无关输入和评分函数的实现
 * 
 * 该文件实现了AutoDock Vina中用于计算分子构象无关特征的功能，
 * 包括扭转角数量、重原子数量、疏水原子数量等特征的计算，
 * 以及基于这些特征的多种评分函数。

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

#include "conf_independent.h"
#include "model.h"

/**
 * @brief 将构象无关输入转换为浮点数向量
 * @return flv 包含所有特征值的浮点数向量
 */
conf_independent_inputs::operator flv() const {
    flv tmp;
    tmp.push_back(num_tors);
    tmp.push_back(num_rotors);
    tmp.push_back(num_heavy_atoms);
    tmp.push_back(num_hydrophobic_atoms);
    tmp.push_back(ligand_max_num_h_bonds);
    tmp.push_back(num_ligands);
    tmp.push_back(ligand_lengths_sum);
    return tmp;
}

/**
 * @brief 计算原子连接的重原子数量
 * @param m 分子模型
 * @param i 原子索引
 * @return unsigned 连接的重原子数量（不包括氢原子）
 * @note 重原子指除氢原子外的所有原子
 */
unsigned conf_independent_inputs::num_bonded_heavy_atoms(const model& m, const atom_index& i) const { 
    unsigned acc = 0;
    const std::vector<bond>& bonds = m.get_atom(i).bonds;

    VINA_FOR_IN(j, bonds) {
        const bond& b = bonds[j];
        const atom& a = m.get_atom(b.connected_atom_index);
        if (!a.is_hydrogen())  // 只计算非氢原子
            ++acc;
    }

    return acc;
}

/**
 * @brief 计算原子的可旋转键数量
 * @param m 分子模型
 * @param i 原子索引
 * @return unsigned 连接到该原子的可旋转键数量
 * @note 只计算连接到重配体原子的可旋转键，排除甲基等末端基团
 */
unsigned conf_independent_inputs::atom_rotors(const model& m, const atom_index& i) const { 
    unsigned acc = 0;
    const std::vector<bond>& bonds = m.get_atom(i).bonds;

    VINA_FOR_IN(j, bonds) {
        const bond& b = bonds[j];
        const atom& a = m.get_atom(b.connected_atom_index);
        // 条件：可旋转 && 非氢原子 && 连接的重原子数>1（排除CH3等末端基团）
        if (b.rotatable && !a.is_hydrogen() && num_bonded_heavy_atoms(m, b.connected_atom_index) > 1) { 
            ++acc;
        }
    }

    return acc;
}

/**
 * @brief 默认构造函数
 * @note 将所有特征值初始化为0
 */
conf_independent_inputs::conf_independent_inputs() : 
    num_tors(0), num_rotors(0), num_heavy_atoms(0), 
    num_hydrophobic_atoms(0), ligand_max_num_h_bonds(0), num_ligands(0), 
    ligand_lengths_sum(0) { }

/**
 * @brief 从分子模型构造构象无关输入
 * @param m 分子模型
 * @note 计算所有配体的构象无关特征，包括扭转角、旋转键、重原子等
 */
conf_independent_inputs::conf_independent_inputs(const model& m) {
    // 初始化所有计数器
    torsdof = 0;
    num_tors = 0;
    num_rotors = 0;
    num_heavy_atoms = 0;
    num_hydrophobic_atoms = 0;
    ligand_max_num_h_bonds = 0;
    num_ligands = m.num_ligands();
    ligand_lengths_sum = 0;

    // 遍历所有配体
    VINA_FOR(i, num_ligands) {
        const ligand& lig = m.get_ligand(i);
        ligand_lengths_sum += m.ligand_length(i);  // 累加配体长度
        torsdof += lig.degrees_of_freedom;         // 累加自由度

        // 遍历配体中的所有原子
        VINA_RANGE(j, lig.begin, lig.end) {
            const atom& a = m.get_atom(j);

            if (a.el != EL_TYPE_H) {  // 只处理非氢原子
                unsigned ar = atom_rotors(m, atom_index(j, false));
                num_tors += 0.5 * ar;  // 扭转角数量（每个键贡献0.5）
                
                // 旋转键数量计算
                if (ar > 2) num_rotors += 0.5;
                else num_rotors += 0.5 * ar;

                // 疏水原子计数
                if (xs_is_hydrophobic(a.xs))
                    ++num_hydrophobic_atoms;

                // 氢键供体/受体计数
                if (xs_is_acceptor(a.xs) || xs_is_donor(a.xs))
                    ++ligand_max_num_h_bonds;

                ++num_heavy_atoms;  // 重原子计数
            }
        }
    }
}

/**
 * @brief 获取特征名称列表
 * @return std::vector<std::string> 特征名称的字符串向量
 * @note 顺序与operator flv()中的顺序一致
 */
std::vector<std::string> conf_independent_inputs::get_names() const { 
    std::vector<std::string> tmp;
    tmp.push_back("num_tors");
    tmp.push_back("num_rotors");
    tmp.push_back("num_heavy_atoms");
    tmp.push_back("num_hydrophobic_atoms");
    tmp.push_back("ligand_max_num_h_bonds");
    tmp.push_back("num_ligands");
    tmp.push_back("ligand_lengths_sum");
    VINA_CHECK(static_cast<flv>(*this).size() == tmp.size()); // 确保大小一致
    return tmp;
}

/**
 * @brief 从迭代器读取下一个浮点数值
 * @param i 浮点数向量的常量迭代器引用
 * @return fl 当前迭代器指向的值
 * @note 读取后迭代器会自动前进
 */
inline fl read_iterator(flv::const_iterator& i) {
    fl x = *i; 
    ++i;
    return x;
}

/**
 * @brief 平滑除法函数，避免除零错误
 * @param x 被除数
 * @param y 除数
 * @return fl 除法结果，对极值情况进行处理
 * @note 当除数接近零时返回适当的极值，避免数值不稳定
 */
fl conf_smooth_div(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); 
    return x / y;
}

// ================ Vina评分函数实现 ================

/**
 * @brief 扭转角平方项评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加基于扭转角数量平方的惩罚项，权重范围[-0.1, 0.1]
 */
fl num_tors_sqr::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.1 * read_iterator(i); // 权重范围[-0.1, 0.1]
    return x + weight * sqr(fl(in.num_tors)) / 5;
}

/**
 * @brief 扭转角平方根项评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加基于扭转角数量平方根的项，相比平方项增长更缓慢
 */
fl num_tors_sqrt::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.1 * read_iterator(i); // 权重范围[-0.1, 0.1]
    return x + weight * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
}

/**
 * @brief 扭转角除法项评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 使用平滑除法，分母包含扭转角数量，权重范围[0, 0.2]
 */
fl num_tors_div::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.1 * (read_iterator(i) + 1); // 权重范围[0, 0.2]
    return conf_smooth_div(x, 1 + weight * in.num_tors / 5.0);
}

/**
 * @brief 配体长度评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加与配体总长度成比例的项
 */
fl ligand_length::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = read_iterator(i);
    return x + weight * in.ligand_lengths_sum;
}

/**
 * @brief 配体数量评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加与配体数量成比例的项，权重范围[-1, 1]
 */
fl num_ligands::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 1 * read_iterator(i); // 权重范围[-1, 1]
    return x + weight * in.num_ligands;
}

/**
 * @brief 重原子数量除法项评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 使用平滑除法，分母包含重原子数量
 */
fl num_heavy_atoms_div::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.05 * read_iterator(i); 
    return conf_smooth_div(x, 1 + weight * in.num_heavy_atoms); 
}

/**
 * @brief 重原子数量评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加与重原子数量成比例的项
 */
fl num_heavy_atoms::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.05 * read_iterator(i); 
    return x + weight * in.num_heavy_atoms;
}

/**
 * @brief 疏水原子数量评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 添加与疏水原子数量成比例的项，体现疏水相互作用
 */
fl num_hydrophobic_atoms::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = 0.05 * read_iterator(i); 
    return x + weight * in.num_hydrophobic_atoms;
}

/**
 * @brief AutoDock 4.2扭转角加法项评分函数
 * @param in 构象无关输入
 * @param x 当前评分值
 * @param i 权重参数迭代器
 * @return fl 更新后的评分值
 * @note 使用AutoDock 4.2中的TORSDOF值，与Vina的扭转角计算方法不同
 */
fl ad4_tors_add::eval(const conf_independent_inputs& in, fl x, flv::const_iterator& i) {
    fl weight = read_iterator(i);
    return x + weight * in.torsdof;
}
