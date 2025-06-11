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

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <boost/optional.hpp> // for context

#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "igrid.h"
#include "grid_dim.h"

/**
 * @brief 原子间相互作用对结构
 * 
 * 存储两个原子之间的相互作用信息，用于计算分子内和分子间的能量贡献
 */
struct interacting_pair {
    sz type_pair_index; ///< 原子类型对的索引，用于查找相互作用参数
    sz a;              ///< 第一个原子的索引
    sz b;              ///< 第二个原子的索引

    /**
     * @brief 构造函数
     * @param type_pair_index_ 原子类型对索引
     * @param a_ 第一个原子索引
     * @param b_ 第二个原子索引
     */
    interacting_pair(sz type_pair_index_, sz a_, sz b_) : type_pair_index(type_pair_index_), a(a_), b(b_) {}
};

typedef std::vector<interacting_pair> interacting_pairs; ///< 相互作用对的容器类型

typedef std::pair<std::string, boost::optional<sz> > parsed_line; ///< 解析后的PDBQT文件行
typedef std::vector<parsed_line> context; ///< PDBQT文件上下文，保存原始文件信息

/**
 * @brief 配体结构
 * 
 * 继承自flexible_body和atom_range，表示分子对接中的小分子配体
 * 包含自由度信息、内部相互作用对和原始文件上下文
 */
struct ligand : public flexible_body, atom_range {
    unsigned degrees_of_freedom; ///< 自由度数量，可能与表观可旋转键数不同（因为禁用的扭转）
    interacting_pairs pairs;     ///< 配体内部的相互作用对
    context cont;               ///< 配体的PDBQT文件上下文
    
    /**
     * @brief 构造函数
     * @param f 柔性体结构
     * @param degrees_of_freedom_ 自由度数量
     */
    ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
    
    /**
     * @brief 设置原子范围
     */
    void set_range();
};

/**
 * @brief 受体残基结构
 * 
 * 继承自main_branch，表示受体蛋白质的柔性残基
 */
struct residue : public main_branch {
    /**
     * @brief 构造函数
     * @param m 主分支结构
     */
    residue(const main_branch& m) : main_branch(m) {}
};

/**
 * @brief 距离类型枚举
 * 
 * 定义原子间距离的变化特性，用于确定相互作用类型
 */
enum distance_type {
    DISTANCE_FIXED,    ///< 固定距离（如共价键）
    DISTANCE_ROTOR,    ///< 转子距离（可旋转键）
    DISTANCE_VARIABLE  ///< 可变距离（非键相互作用）
};

typedef strictly_triangular_matrix<distance_type> distance_type_matrix; ///< 距离类型矩阵

struct cache;               ///< 缓存结构
struct szv_grid;            ///< 网格结构
struct pdbqt_initializer;   ///< PDBQT初始化器（仅在parse_pdbqt.cpp中声明）
struct precalculate_byatom; ///< 按原子预计算结构

/**
 * @brief 计算相互作用对的能量
 * @param p 按原子预计算的相互作用参数
 * @param v 当前能量值
 * @param pairs 相互作用对列表
 * @param coords 原子坐标
 * @param with_max_cutoff 是否应用最大截断距离
 * @return 计算得到的能量值
 */
fl eval_interacting_pairs(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, const bool with_max_cutoff=false);

/**
 * @brief 计算相互作用对的能量和导数（力）
 * @param p 按原子预计算的相互作用参数
 * @param v 当前能量值
 * @param pairs 相互作用对列表
 * @param coords 原子坐标
 * @param forces 输出的原子受力
 * @param with_max_cutoff 是否应用最大截断距离
 * @return 计算得到的能量值
 */
fl eval_interacting_pairs_deriv(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces, const bool with_max_cutoff=false);

/**
 * @brief 分子模型核心类
 * 
 * 表示完整的分子系统，包含配体、受体柔性部分、原子坐标、相互作用等信息
 * 是AutoDock Vina分子对接算法的核心数据结构
 */
struct model {
private:
    // 友元类声明 - 允许这些类访问私有成员
    friend struct cache;
    friend struct non_cache;
    friend struct ad4cache;
    friend struct szv_grid;
    friend struct appender_info;
    friend struct pdbqt_initializer;

    // 私有文件输出方法
    void write_context(const context &c, std::ostringstream& out) const; ///< 写入上下文到字符串流
    void write_context(const context& c, ofile& out) const; ///< 写入上下文到文件
    
    /**
     * @brief 写入上下文到文件（带备注）
     * @param c 上下文
     * @param out 输出文件
     * @param remark 备注
     */
    void write_context(const context& c, ofile& out, const std::string& remark) const {
        out << remark;
    }
    
    /**
     * @brief 写入上下文到文件路径
     * @param c 上下文
     * @param name 文件路径
     */
    void write_context(const context& c, const path& name) const {
        ofile out(name);
        write_context(c, out);
    }
    
    /**
     * @brief 写入上下文到文件路径（带备注）
     * @param c 上下文
     * @param name 文件路径
     * @param remark 备注
     */
    void write_context(const context& c, const path& name, const std::string& remark) const {
        ofile out(name);
        write_context(c, out, remark);
    }
    
    /**
     * @brief 计算非对称RMSD下界
     * @param x 第一个模型
     * @param y 第二个模型
     * @return RMSD下界
     * @note 实际上是静态方法
     */
    fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const;
    
    // 私有辅助方法
    atom_index sz_to_atom_index(sz i) const; ///< 大小索引转换为原子索引（grid_atoms, atoms）
    bool bonded_to_HD(const atom& a) const; ///< 判断是否与氢原子键合
    bool bonded_to_heteroatom(const atom& a) const; ///< 判断是否与杂原子键合
    void bonded_to(sz a, sz n, szv& out) const; ///< 查找键合原子
    szv bonded_to(sz a, sz n) const; ///< 查找键合原子（返回版本）
    bool is_closure_clash(sz i, sz j) const; ///< 判断是否为闭环碰撞
    bool is_unmatched_closure_dummy(sz i, sz j) const; ///< 判断是否为不匹配的闭环虚原子
    bool is_glue_pair(sz i, sz j) const; ///< 判断是否为胶合对

    /**
     * @brief 基于相对迁移率分配键
     * @param mobility 迁移率矩阵
     * @note 基于相对迁移率、距离和共价长度分配键
     */
    void assign_bonds(const distance_type_matrix& mobility);
    
    void assign_types(); ///< 分配原子类型
    
    /**
     * @brief 初始化相互作用对
     * @param mobility 迁移率矩阵
     */
    void initialize_pairs(const distance_type_matrix& mobility);
    
    /**
     * @brief 初始化模型
     * @param mobility 迁移率矩阵
     */
    void initialize(const distance_type_matrix& mobility);
    
    /**
     * @brief 计算碰撞惩罚辅助函数
     * @param pairs 相互作用对
     * @return 碰撞惩罚值
     */
    fl clash_penalty_aux(const interacting_pairs& pairs) const;

    // 核心数据成员
    vecv coords;        ///< 原子坐标向量
    vecv minus_forces;  ///< 负力向量

    atomv grid_atoms;   ///< 网格原子
    atomv atoms;        ///< 可移动和不可弯曲原子
    vector_mutable<ligand> ligands;  ///< 配体列表
    vector_mutable<residue> flex;    ///< 柔性残基列表
    context flex_context;            ///< 柔性部分的PDBQT上下文
    
    interacting_pairs other_pairs;   ///< 分子内相互作用：flex_i - flex_j 和 flex_i - flex_i
    interacting_pairs inter_pairs;   ///< 分子间相互作用：ligand - flex 和 ligand_i - ligand_j
    interacting_pairs glue_pairs;    ///< 分子内相互作用：glue_i - glue_i

    sz m_num_movable_atoms;          ///< 可移动原子数量
    atom_type::t m_atom_typing_used; ///< 使用的原子类型系统
public:
    /**
     * @brief 默认构造函数
     * @note 初始化可移动原子数为0，原子类型为XS
     */
    model() : m_num_movable_atoms(0), m_atom_typing_used(atom_type::XS) {}
    
    /**
     * @brief 指定原子类型的构造函数
     * @param atype 使用的原子类型
     */
    model(atom_type::t atype) : m_num_movable_atoms(0), m_atom_typing_used(atype) {}

    // 数据访问方法
    atomv get_atoms() const { return atoms; }               ///< 获取所有原子（用于precalculate_byatom）
    atom get_atom(sz i) const { return atoms[i]; }          ///< 获取指定原子
    ligand get_ligand(sz i) const { return ligands[i]; }    ///< 获取指定配体
    vec get_coords(sz i) const { return coords[i]; }        ///< 获取指定原子坐标
    interacting_pairs get_other_pairs() const { return other_pairs; }   ///< 获取其他相互作用对
    interacting_pairs get_inter_pairs() const { return inter_pairs; }   ///< 获取分子间相互作用对
    interacting_pairs get_glue_pairs() const { return glue_pairs; }     ///< 获取胶合相互作用对

    /**
     * @brief 合并另一个模型到当前模型
     * @param m 要合并的模型
     */
    void append(const model& m);
    
    /**
     * @brief 获取使用的原子类型
     * @return 原子类型
     */
    atom_type::t atom_typing_used() const { return m_atom_typing_used; }

    // 查询方法
    bool is_atom_in_ligand(sz a) const; ///< 判断原子是否在配体中
    bool is_movable_atom(sz a) const; ///< 判断原子是否可移动
    std::vector<double> center() const; ///< 计算几何中心
    sz find_ligand(sz a) const; ///< 查找原子所属的配体
    sz num_atoms() const { return atoms.size(); } ///< 原子总数
    sz num_movable_atoms() const { return m_num_movable_atoms; } ///< 可移动原子数
    sz num_internal_pairs() const; ///< 内部相互作用对数
    sz num_other_pairs() const { return other_pairs.size(); } ///< 其他相互作用对数
    sz num_ligands() const { return ligands.size(); } ///< 配体数量
    sz num_flex() const { return flex.size(); } ///< 柔性残基数量
    sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; } ///< 配体自由度
    sz ligand_longest_branch(sz ligand_number) const; ///< 配体最长分支
    sz ligand_length(sz ligand_number) const; ///< 配体长度

    /**
     * @brief 获取可移动原子的类型
     * @param atom_typing_used_ 使用的原子类型系统
     * @return 原子类型向量
     */
    szv get_movable_atom_types(atom_type::t atom_typing_used_) const;
    
    vecv get_ligand_coords() const; ///< 获取配体坐标
    std::vector<double> get_ligand_coords(); ///< 获取配体坐标（双精度版本）
    vecv get_heavy_atom_movable_coords() const; ///< 获取重原子可移动坐标
    conf_size get_size() const; ///< 获取构象大小
    conf get_initial_conf() const; ///< 获取初始构象（扭转=0，取向=单位，配体位置=当前）

    // 文件输出方法
    /**
     * @brief 写入柔性部分
     * @param name 文件路径
     * @param remark 备注信息
     */
    void write_flex(const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }
    
    /**
     * @brief 写入指定配体
     * @param ligand_number 配体编号
     * @param name 文件路径
     * @param remark 备注信息
     */
    void write_ligand(sz ligand_number, const path& name, const std::string& remark) const { 
        VINA_CHECK(ligand_number < ligands.size()); 
        write_context(ligands[ligand_number].cont, name, remark); 
    }
    
    /**
     * @brief 写入结构到输出文件
     * @param out 输出文件流
     */
    void write_structure(ofile& out) const {
        VINA_FOR_IN(i, ligands)
            write_context(ligands[i].cont, out);
        if(num_flex() > 0) // 否则备注会被徒劳写入
            write_context(flex_context, out);
    }
    
    /**
     * @brief 写入结构到输出文件（带备注）
     * @param out 输出文件流
     * @param remark 备注信息
     */
    void write_structure(ofile& out, const std::string& remark) const {
        out << remark;
        write_structure(out);
    }
    
    /**
     * @brief 写入结构到输出文件（带多行备注）
     * @param out 输出文件流
     * @param remarks 备注信息列表
     */
    void write_structure(ofile& out, std::vector<std::string>& remarks) const {
        VINA_FOR_IN(i, remarks) out << remarks[i];
        write_structure(out);
    }

    /**
     * @brief 写入结构到文件
     * @param name 文件路径
     */
    void write_structure(const path& name) const { ofile out(name); write_structure(out); }
    
    /**
     * @brief 写入模型到输出文件
     * @param out 输出文件流
     * @param model_number 模型编号
     * @param remark 备注信息
     */
    void write_model(ofile& out, sz model_number, const std::string& remark) const {
        out << "MODEL " << model_number << '\n';
        write_structure(out, remark);
        out << "ENDMDL\n";
    }
    
    /**
     * @brief 写入模型并返回字符串
     * @param model_number 模型编号
     * @param remark 备注信息
     * @return 模型的字符串表示
     */
    std::string write_model(sz model_number, const std::string &remark);

    /**
     * @brief 设置模型构象
     * @param c 构象参数
     */
    void set(const conf& c);

    /**
     * @brief 计算配体的回转半径
     * @param ligand_number 配体编号
     * @return 回转半径
     * @note 使用当前坐标计算
     */
    fl gyration_radius(sz ligand_number) const;

    // 原子访问方法
    const atom_base& movable_atom(sz i) const { assert(i < m_num_movable_atoms); return atoms[i]; } ///< 获取可移动原子
    const vec& movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; } ///< 获取可移动原子坐标

    const vec& atom_coords(const atom_index& i) const; ///< 获取原子坐标
    fl distance_sqr_between(const atom_index& a, const atom_index& b) const; ///< 计算原子间距离平方
    
    /**
     * @brief 判断两原子间是否存在中间原子
     * @param mobility 迁移率矩阵
     * @param a 第一个原子
     * @param b 第二个原子
     * @param relevant_atoms 相关原子列表
     * @return 是否存在中间原子
     * @note 存在一个原子同时距离a和b更近，且相对于它们不可移动
     */
    bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const;

    /**
     * @brief 获取两原子间的距离类型
     * @param mobility 迁移率矩阵
     * @param i 第一个原子索引
     * @param j 第二个原子索引
     * @return 距离类型
     */
    distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

    // 能量计算方法
    fl evalo(const precalculate_byatom& p, const vec& v) const; ///< 计算外部能量
    fl evali(const precalculate_byatom& p, const vec& v) const; ///< 计算内部能量
    fl eval_inter(const precalculate_byatom& p, const vec& v) const; ///< 计算分子间能量
    
    /**
     * @brief 计算能量导数
     * @param p 预计算参数
     * @param ig 交互网格
     * @param v 位置向量
     * @param g 梯度变化
     * @return 能量值
     */
    fl eval_deriv(const precalculate_byatom& p, const igrid& ig, const vec& v, change& g);
    
    /**
     * @brief 计算分子内能量
     * @param p 预计算参数
     * @param ig 交互网格
     * @param v 位置向量
     * @return 分子内能量
     */
    fl eval_intramolecular(const precalculate_byatom& p, const igrid& ig, const vec& v);

    // RMSD计算方法
    fl rmsd_lower_bound(const model& m) const; ///< RMSD下界（使用coords）
    fl rmsd_upper_bound(const model& m) const; ///< RMSD上界（使用coords）
    fl rmsd_ligands_upper_bound(const model& m) const; ///< 配体RMSD上界（使用coords）

    // 验证和调试方法
    void verify_bond_lengths() const; ///< 验证键长
    void about() const; ///< 显示模型信息
    void check_ligand_internal_pairs() const; ///< 检查配体内部相互作用对
    void show_pairs() const; ///< 显示相互作用对
    void show_atoms() const; ///< 显示原子信息
    void show_forces() const; ///< 显示力信息
    
    /**
     * @brief 打印调试信息
     * @param show_coords 是否显示坐标
     * @param show_internal 是否显示内部信息
     * @param show_atoms 是否显示原子信息
     * @param show_grid 是否显示网格信息
     * @param show_about 是否显示关于信息
     * @note FIXME 后续可能需要移除
     */
    void print_stuff(bool show_coords=true, bool show_internal=true, bool show_atoms=true, bool show_grid=true, bool show_about=true) const;

    /**
     * @brief 计算碰撞惩罚
     * @return 碰撞惩罚值
     */
    fl clash_penalty() const;

    // 原子获取方法
    const atom& get_atom(const atom_index& i) const { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); } ///< 获取原子（常量版本）
    atom& get_atom(const atom_index& i) { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); } ///< 获取原子（可修改版本）
};

#endif
