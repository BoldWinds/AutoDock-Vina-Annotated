/**
 * @file parse_pdbqt.cpp
 * @brief PDBQT解析和模型构建实现

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

#include <string>
#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <exception>
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
//#include <openbabel/mol.h>
//#include <openbabel/obconversion.h>
#include "model.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "utils.h"
#include "parse_pdbqt.h"
#include "parse_error.h"


/**
 * @struct parsed_atom
 * @brief 解析后的原子结构体，保存原子序号；继承自基础atom类，添加原子类型、电荷与坐标信息
 */
struct parsed_atom : public atom {
    unsigned number; 
    
    /**
     * @brief 构造函数
     * @param ad_ AutoDock原子类型ID
     * @param charge_ 部分电荷
     * @param coords_ 三维坐标
     * @param number_ PDBQT文件中的原子序号
     * 
     * @note 该构造函数初始化原子的所有基本属性，为后续的
     * 分子结构构建和键连关系建立提供基础数据
     */
    parsed_atom(sz ad_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
        ad = ad_;
        charge = charge_;
        coords = coords_;
    }
};

/**
 * @brief 向上下文中添加解析行信息
 * @param c 上下文对象的引用
 * @param str 待添加的字符串行
 */
void add_context(context& c, std::string& str) {
    c.push_back(parsed_line(str, boost::optional<sz>()));
}

/**
 * @brief 从字符串的指定范围内提取内容并去除首尾空白字符
 * @param str 源字符串
 * @param i 起始位置（基于1的索引）
 * @param j 结束位置（基于1的索引，包含该位置）
 * @return std::string 去除空白字符后的子字符串
 * 
 * @note 似乎有bug，j>=str.size()一直成立
 */
std::string omit_whitespace(const std::string& str, sz i, sz j) {
    if(i < 1) i = 1;
    if(j < i-1) j = i-1; // i >= 1
    if(j < str.size()) j = str.size();

    // 去除前导空白字符（从左向右扫描）
    while(i <= j && std::isspace(str[i-1]))
        ++i;

    // 去除尾随空白字符（从右向左扫描）
    while(i <= j && std::isspace(str[j-1]))
        --j;

    // 安全性检查：确保计算出的索引有效
    VINA_CHECK(i-1 < str.size());      // 起始索引有效
    VINA_CHECK(j-i+1 <= str.size());   // 子字符串长度有效

    // 提取并返回处理后的子字符串（转换为基于0的索引）
    return str.substr(i-1, j-i+1);
}

/**
 * @brief 安全地从字符串的指定位置提取并转换为指定类型的数值
 * @tparam T 目标转换类型（如unsigned、fl等）
 * @param str 源字符串
 * @param i 起始位置（基于1的索引）
 * @param j 结束位置（基于1的索引，包含该位置）
 * @param dest_nature 目标数据的性质描述，用于错误信息
 * @return T 转换后的目标类型数值
 * @throws pdbqt_parse_error 当字符串过短或转换失败时抛出异常
 */
template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
    // 验证索引参数的有效性
    VINA_CHECK(i >= 1);
    VINA_CHECK(i <= j+1);

    if(j > str.size()) 
        throw pdbqt_parse_error("This line is too short.", str);

    // 跳过前导空白字符（PDBQT格式中字段可能有前导空格）
    while(i <= j && std::isspace(str[i-1]))
        ++i;

    // 提取子字符串（转换为基于0的索引）
    const std::string substr = str.substr(i-1, j-i+1);
    
    // 尝试类型转换
    try {
        return boost::lexical_cast<T>(substr);
    }
    catch(...) {
        throw pdbqt_parse_error(dest_nature + std::string(" \"") + substr + "\" is not valid.", str);
    }
}

/**
 * @brief 解析PDBQT格式的原子行字符串
 * @param str PDBQT文件中的原子行字符串（ATOM或HETATM行）
 * @return parsed_atom 解析后的原子对象
 * @throws pdbqt_parse_error 当原子信息格式错误或原子类型无效时抛出异常
 * 
 * @note PDBQT原子字段举例：
 * ATOM      1  C1  LIG     0       1.207   0.000   0.000  -0.001 C
 * - 位置7-11：原子序号
 * - 位置31-38：X坐标
 * - 位置39-46：Y坐标  
 * - 位置47-54：Z坐标
 * - 位置69-76：部分电荷（可选）
 * - 位置78-79：AutoDock原子类型名称
 */
parsed_atom parse_pdbqt_atom_string(const std::string& str) {
    // 解析原子序号（位置7-11）
    unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "Atom number");
    
    // 解析三维坐标（位置31-54）
    vec coords(checked_convert_substring<fl>(str, 31, 38, "Coordinate"),
               checked_convert_substring<fl>(str, 39, 46, "Coordinate"),
               checked_convert_substring<fl>(str, 47, 54, "Coordinate"));
    
    // 解析部分电荷（位置69-76，可选字段）
    fl charge = 0;
    if(!substring_is_blank(str, 69, 76))  // 检查电荷字段是否为空
        charge = checked_convert_substring<fl>(str, 69, 76, "Charge");
    
    // 解析原子类型名称（位置78-79）并去除空白字符
    std::string name = omit_whitespace(str, 78, 79);
    
    // 将原子类型名称转换为AutoDock内部类型ID
    sz ad = string_to_ad_type(name);
    
    // 创建解析后的原子对象
    parsed_atom tmp(ad, charge, coords, number);

    // 特殊处理：如果是非AutoDock标准的金属原子，设置为金属类型
    if(is_non_ad_metal_name(name))
        tmp.xs = XS_TYPE_Met_D;  // 设置为金属供体类型
    
    // 验证原子类型的有效性
    if(tmp.acceptable_type()) 
        return tmp;  // 原子类型有效，返回解析结果
    else
        // 原子类型无效，抛出详细错误信息
        throw pdbqt_parse_error("Atom type " + name + " is not a valid AutoDock type (atom types are case-sensitive).", str);
}

/**
 * @struct atom_reference
 * @brief 原子引用结构体，用于在解析过程中跟踪原子的位置和类型
 * 
 * 在PDBQT解析过程中，需要区分可移动原子和不可移动原子，此结构体
 * 提供了统一的引用方式来标识原子在不同数组中的位置。
 */
struct atom_reference {
    sz index;    ///< 原子在相应数组中的索引位置
    bool inflex; ///< 是否为不可移动原子（true=不可移动，false=可移动）
    
    /**
     * @brief 构造函数
     * @param index_ 原子索引
     * @param inflex_ 是否为不可移动原子
     */
    atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

/**
 * @struct movable_atom
 * @brief 可移动原子结构体，继承自基础atom类并添加相对坐标信息
 * 
 * 在分子对接过程中，可移动原子需要保存相对于其所属刚体的相对坐标，
 * 以便在旋转和平移变换时正确计算新的原子位置。
 */
struct movable_atom : public atom {
    vec relative_coords; ///< 相对于所属刚体原点的相对坐标
    
    /**
     * @brief 构造函数
     * @param a 基础原子对象
     * @param relative_coords_ 相对坐标
     */
    movable_atom(const atom& a, const vec& relative_coords_) : atom(a) {
        relative_coords = relative_coords_;
    }
};

/**
 * @struct rigid
 * @brief 刚性结构体，包含不可移动的原子集合
 * 
 * 代表受体的刚性部分，这些原子在对接过程中位置固定不变
 */
struct rigid {
    atomv atoms; ///< 刚性原子向量
};

typedef std::vector<movable_atom> mav; ///< 可移动原子向量的类型别名

/**
 * @struct non_rigid_parsed
 * @brief 非刚性解析结果结构体，包含所有解析后的可移动分子组件
 * 
 * 该结构体是PDBQT解析的核心数据容器，包含了配体、柔性残基以及
 * 原子间的键连关系矩阵。这些信息将用于构建最终的分子模型。
 */
struct non_rigid_parsed {
    vector_mutable<ligand> ligands;  ///< 配体向量
    vector_mutable<residue> flex;    ///< 柔性残基向量

    mav atoms;     ///< 可移动原子向量
    atomv inflex;  ///< 不可移动原子向量

    distance_type_matrix atoms_atoms_bonds;   ///< 可移动原子间的键连关系矩阵
    matrix<distance_type> atoms_inflex_bonds; ///< 可移动原子与不可移动原子间的键连关系矩阵
    distance_type_matrix inflex_inflex_bonds; ///< 不可移动原子间的键连关系矩阵

    /**
     * @brief 获取完整的移动性矩阵
     * @return distance_type_matrix 包含所有原子间移动性约束的矩阵
     * 
     * @note 将三个分块矩阵合并为一个完整矩阵，用于后续的力场计算
     */
    distance_type_matrix mobility_matrix() const {
        distance_type_matrix tmp(atoms_atoms_bonds);
        tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
        return tmp;
    }
};

/**
 * @struct parsing_struct
 * @brief 解析结构体，用于在PDBQT解析过程中维护层次化的分子结构
 * 
 * PDBQT文件使用ROOT/BRANCH结构描述分子的层次化组织，此结构体
 * 递归地表示这种树形结构，支持任意深度的分支嵌套。
 */
struct parsing_struct {
    /**
     * @struct node_t
     * @brief 解析树节点模板结构体，表示分子树中的一个节点
     * @tparam T 通常为parsing_struct类型，支持递归定义
     * 
     * 每个节点包含一个原子和可能的子分支，构成了分子的完整树形结构。
     */
    template<typename T> // T == parsing_struct
    struct node_t {
        sz context_index;    ///< 在上下文数组中的索引位置
        parsed_atom a;       ///< 节点包含的原子
        std::vector<T> ps;   ///< 子分支向量
        
        /**
         * @brief 构造函数
         * @param a_ 原子对象
         * @param context_index_ 上下文索引
         */
        node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

        /**
         * @brief 将该节点的所有子分支的axis_begin（轴开始）设置为当前节点的原子引用；
         *        并把当前原子插入到nr的不可移动原子集合中
         * @param nr 非刚性解析结果引用
         */
        void insert_inflex(non_rigid_parsed& nr) {
            VINA_FOR_IN(i, ps)
                ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
            nr.inflex.push_back(a);
        }
        
        /**
         * @brief 递归的设置所有子分支的轴结束
         * @param nr 非刚性解析结果引用
         */
        void insert_immobiles_inflex(non_rigid_parsed& nr) {
            VINA_FOR_IN(i, ps)
                ps[i].insert_immobile_inflex(nr);
        }

        /**
         * @brief 将原子插入到可移动原子集合中(用相对位置)，并给存在的子分支设置axis_begin
         * @param nr 非刚性解析结果引用
         * @param c 上下文引用
         * @param frame_origin 参考坐标系原点
         * 
         * @details 对于配体，frame_origin就是根的第一个原子的坐标；
         *          对于柔性残基，frame_origin是该分支的immobile_atom坐标
         */
        void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
            VINA_FOR_IN(i, ps)
                ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
            vec relative_coords; relative_coords = a.coords - frame_origin;
            c[context_index].second = nr.atoms.size();
            nr.atoms.push_back(movable_atom(a, relative_coords));
        }
        
        /**
         * @brief 递归插入所有不可移动的子原子
         * @param nr 非刚性解析结果引用
         * @param c 上下文引用
         * @param frame_origin 参考坐标系原点
         */
        void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
            VINA_FOR_IN(i, ps)
                ps[i].insert_immobile(nr, c, frame_origin);
        }
    };

    typedef node_t<parsing_struct> node; ///< 节点类型别名
    
    boost::optional<sz> immobile_atom;                ///< 不可移动原子的索引（本分支与父分支相连的原子）
    boost::optional<atom_reference> axis_begin;       ///< 轴起始原子引用（由父节点设置）
    boost::optional<atom_reference> axis_end;         ///< 轴结束原子引用
    std::vector<node> atoms;                          ///< 节点向量

    /**
     * @brief 添加原子到解析结构中
     * @param a 解析后的原子
     * @param c 上下文
     */
    void add(const parsed_atom& a, const context& c) { 
        VINA_CHECK(c.size() > 0);
        atoms.push_back(node(a, c.size()-1)); 
    }
    
    /**
     * @brief 获取不可移动原子的坐标
     * @return const vec& 不可移动原子的坐标引用
     */
    const vec& immobile_atom_coords() const {
        VINA_CHECK(immobile_atom);
        VINA_CHECK(immobile_atom.get() < atoms.size());
        return atoms[immobile_atom.get()].a.coords;
    }
    
    /**
     * @brief 设置axis_end的原子引用并将immobile_atom插入不可移动原子集合
     * @param nr 非刚性解析结果引用
     */
    void insert_immobile_inflex(non_rigid_parsed& nr) {
        if(!atoms.empty()) {
            VINA_CHECK(immobile_atom);
            VINA_CHECK(immobile_atom.get() < atoms.size());
            axis_end = atom_reference(nr.inflex.size(), true);
            atoms[immobile_atom.get()].insert_inflex(nr);
        }
    }

    /**
     * @brief 若该分支有原子，则设置axis_end的原子引用，并继续遍历添加可移动原子
     * @param nr 非刚性解析结果引用
     * @param c 上下文引用
     * @param frame_origin 参考坐标系原点
     * 
     * @note 如果存在原子，设置轴结束引用并插入可移动原子
     */
    void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
        if(!atoms.empty()) {
            VINA_CHECK(immobile_atom);
            VINA_CHECK(immobile_atom.get() < atoms.size());
            axis_end = atom_reference(nr.atoms.size(), false);
            atoms[immobile_atom.get()].insert(nr, c, frame_origin);
        }
    }

    /**
     * @brief 检查该节点是否只包含一个不可变的原子且该原子不接入分支
     * @return bool 检查结果
     * @note 如果检查结果为true则说明没有必要为该节点建立分支
     */
    bool essentially_empty() const {
        VINA_FOR_IN(i, atoms) {
            if(immobile_atom && immobile_atom.get() != i)
                return false;
            const node& nd = atoms[i];
            if(!nd.ps.empty())
                return false; // FIXME : iffy
        }
        return true;
    }
};

/**
 * @brief 从字符串解析单个无符号整数
 * @param str 源字符串
 * @param start 起始关键字
 * @return unsigned 解析出的无符号整数
 * @throws pdbqt_parse_error 解析失败时抛出异常
 */
unsigned parse_one_unsigned(const std::string& str, const std::string& start) {
    std::istringstream in_str(str.substr(start.size()));
    int tmp;
    in_str >> tmp;

    if(!in_str || tmp < 0) 
        throw pdbqt_parse_error("Syntax error.", str);

    return unsigned(tmp);
}

/**
 * @brief 从字符串解析两个无符号整数
 * @param str 源字符串
 * @param start 起始关键字
 * @param first 第一个整数的引用
 * @param second 第二个整数的引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 */
void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned& first, unsigned& second) {
    std::istringstream in_str(str.substr(start.size()));
    int tmp1, tmp2;
    in_str >> tmp1;
    in_str >> tmp2;

    if(!in_str || tmp1 < 0 || tmp2 < 0) 
        throw pdbqt_parse_error("Syntax error.", str);

    first = unsigned(tmp1);
    second = unsigned(tmp2);
}

/**
 * @brief 解析刚性受体的PDBQT文件
 * @param name 文件路径
 * @param r 刚性结构引用，用于存储解析结果
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 刚性受体只包含ATOM/HETATM行，不支持BRANCH结构
 */
void parse_pdbqt_rigid(const path& name, rigid& r) {
    ifile in(name);
    std::string str;

    while(std::getline(in, str)) {
        if(str.empty()) {} // 忽略空行
        else if(starts_with(str, "TER")) {} // 忽略TER行
        else if(starts_with(str, "END")) {} // 忽略END行
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM"))
            r.atoms.push_back(parse_pdbqt_atom_string(str));
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in rigid receptor. "
                                    "Only one model can be used for the rigid receptor.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in rigid receptor.", str);
    }
}

/**
 * @brief 解析ROOT部分的辅助函数
 * @param in 输入流
 * @param p 解析结构引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 解析ROOT...ENDROOT之间的原子，这些原子属于分子的主干部分
 */
void parse_pdbqt_root_aux(std::istream& in, parsing_struct& p, context& c) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // 忽略空行
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM"))
            p.add(parse_pdbqt_atom_string(str), c);
        else if(starts_with(str, "ENDROOT"))
            return;
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

/**
 * @brief 解析ROOT部分的主函数
 * @param in 输入流
 * @param p 解析结构引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 查找ROOT关键字并调用辅助函数解析其内容
 */
void parse_pdbqt_root(std::istream& in, parsing_struct& p, context& c) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // 忽略空行
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "ROOT")) {
            parse_pdbqt_root_aux(in, p, c);
            break;
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

void parse_pdbqt_branch(std::istream& in, parsing_struct& p, context& c, unsigned from, unsigned to); // 前向声明

/**
 * @brief 解析BRANCH行的辅助函数
 * @param in 输入流
 * @param str BRANCH行字符串
 * @param p 解析结构引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 解析BRANCH行中的两个原子编号，找到对应的父原子并创建子分支
 */
void parse_pdbqt_branch_aux(std::istream& in, const std::string& str, parsing_struct& p, context& c) {
    unsigned first, second;
    parse_two_unsigneds(str, "BRANCH", first, second);
    sz i = 0;

    // 查找起始原子编号对应的原子
    for(; i < p.atoms.size(); ++i)
        if(p.atoms[i].a.number == first) {
            p.atoms[i].ps.push_back(parsing_struct());
            parse_pdbqt_branch(in, p.atoms[i].ps.back(), c, first, second);
            break;
        }

    if(i == p.atoms.size())
        throw pdbqt_parse_error("Atom number " + std::to_string(first) + " is missing in this branch.", str);
}

/**
 * @brief PDBQT解析的主要辅助函数
 * @param in 输入流
 * @param p 解析结构引用
 * @param c 上下文引用
 * @param torsdof 扭转自由度数量（可选）
 * @param residue 是否为残基模式
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 解析ROOT部分后继续解析BRANCH结构，构建完整的分子层次结构
 */
void parse_pdbqt_aux(std::istream& in, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
    parse_pdbqt_root(in, p, c);

    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // 忽略空行
        if(str[0] == '\0') {} // 忽略另一种空行形式（Windows潜在问题）
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, str, p, c);
        else if(!residue && starts_with(str, "TORSDOF")) {
            if(torsdof)
                throw pdbqt_parse_error("TORSDOF keyword can be defined only once.");
            torsdof = parse_one_unsigned(str, "TORSDOF");
        }
        else if(residue && starts_with(str, "END_RES"))
            return; 
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

/**
 * @brief 在非刚性解析结果中添加键连关系
 * @param nr 非刚性解析结果引用
 * @param atm 原子引用（可选）
 * @param r 原子范围
 * 
 * @note 若atm可移动，则添加DISTANCE_FIXED到可移动原子与不可移动原子间的键连关系矩阵atoms_inflex_bonds
 *      否则添加DISTANCE_FIXED到可移动原子与可移动原子间的键连关系矩阵atoms_atoms_bonds
 */
void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
    if(atm)
        VINA_RANGE(i, r.begin, r.end) {
            atom_reference& ar = atm.get();
            if(ar.inflex) 
                nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED;
            else
                nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED;
        }
}

/**
 * @brief 设置旋转轴
 * @param nr 非刚性解析结果引用
 * @param axis_begin 轴起始原子引用
 * @param axis_end 轴结束原子引用
 * 
 * @note 若axis_start和axis_end都不可移动，则设置inflex_inflex_bonds为DISTANCE_ROTOR；
 * 若axis_end可移动，则根据axis_start是否可以移动，向atoms_inflex_bonds或者atoms_atoms_bonds放置DISTANCE_ROTOR
 */
void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
    if(axis_begin && axis_end) {
        atom_reference& r1 = axis_begin.get();
        atom_reference& r2 = axis_end  .get();
        if(r2.inflex) {
            VINA_CHECK(r1.inflex); // 不允许原子-不可移动原子旋转轴
            nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
        }
        else
            if(r1.inflex)
                nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (可移动原子, 不可移动原子)
            else
                nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
    }
}

typedef std::pair<sz, sz> axis_numbers;           ///< 轴编号对的类型别名
typedef boost::optional<axis_numbers> axis_numbers_option; ///< 可选轴编号对的类型别名

/**
 * @brief 更新非刚性解析结果中的距离矩阵
 * @param nr 非刚性解析结果引用
 * 
 * @note 根据当前原子数量调整所有距离矩阵的大小，并设置默认的距离类型
 */
void nr_update_matrixes(non_rigid_parsed& nr) {
    // axis_begin和axis_end索引的原子不能相对于[b.node.begin, b.node.end)范围移动

    nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
    nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // 第一个索引-不可移动原子，第二个索引-可移动原子
    nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

/**
 * @brief 后处理分支的模板函数
 * @tparam B 分支类型（ligand->flexible_body(heterotree<rigid_body>)/residue->main_branch(heterotree<first_segment>)/segment）
 * @param nr 非刚性解析结果引用
 * @param p 目前节点解析结果
 * @param c 上下文引用
 * @param b 分支对象引用；b.node为
 * 
 * @note 这是分子结构构建的核心函数，将解析树转换为最终的分子表示：
 * 1. 插入原子到相应集合中
 * 2. 更新键连矩阵
 * 3. 设置旋转轴
 * 4. 递归处理子分支
 * 
 * @details 配体时b.node为rigid_body；柔性残基时b.node为first_segment；其他b.node为segment
 */
template<typename B> // B == ligand / residue
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
    // 1. 设置分支节点的原子起始索引
    b.node.begin = nr.atoms.size();

    // 2. 插入可移动原子
    VINA_FOR_IN(i, p.atoms) {
        parsing_struct::node& p_node = p.atoms[i];
        if(p.immobile_atom && i == p.immobile_atom.get()) {}
        else p_node.insert(nr, c, b.node.get_origin());
        p_node.insert_immobiles(nr, c, b.node.get_origin());
    }

    // 设置分支节点的原子结束索引
    b.node.end = nr.atoms.size();

    // 3. 更新键连矩阵
    nr_update_matrixes(nr);

    // 4. 添加轴键连关系
    add_bonds(nr, p.axis_begin, b.node);
    add_bonds(nr, p.axis_end  , b.node);
    set_rotor(nr, p.axis_begin, p.axis_end);

    // 5. 设置同一分支内原子间的固定距离约束
    VINA_RANGE(i, b.node.begin, b.node.end)
        VINA_RANGE(j, i+1, b.node.end)
            nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME

    // 6. 递归处理子分支
    VINA_FOR_IN(i, p.atoms) {
        parsing_struct::node& p_node = p.atoms[i];
        VINA_FOR_IN(j, p_node.ps) {
            parsing_struct& ps = p_node.ps[j];
            if(!ps.essentially_empty()) { // 不可移动原子已插入 // FIXME ?!
                // 创建子分支段
                b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch将分配begin和end
                postprocess_branch(nr, ps, c, b.children.back());
            }
        }
    }
    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

/**
 * @brief 后处理配体
 * @param nr 非刚性解析结果引用
 * @param p 解析结构引用
 * @param c 上下文引用
 * @param torsdof 扭转自由度数量
 * 
 * @note 为配体创建flexible_body对象并调用通用的分支后处理函数
 */
void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
    VINA_CHECK(!p.atoms.empty());
    nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch将分配begin和end
    postprocess_branch(nr, p, c, nr.ligands.back());
    nr_update_matrixes(nr); // FIXME ?
}

/**
 * @brief 后处理柔性残基
 * @param nr 非刚性解析结果引用
 * @param p 根的解析结果引用
 * @param c 上下文引用
 * 
 * @note 子分支的解析结果引用在p.atoms中的对应原子节点的ps字段中
 */
void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
    // 把这个残基的所有不可移动原子插入到不可移动集合；并设置axis_begin和axis_end
    VINA_FOR_IN(i, p.atoms) {
        parsing_struct::node& p_node = p.atoms[i];
        p_node.insert_inflex(nr);
        p_node.insert_immobiles_inflex(nr);
    }
    // 处理每一个分支
    VINA_FOR_IN(i, p.atoms) {
        parsing_struct::node& p_node = p.atoms[i];
        VINA_FOR_IN(j, p_node.ps) {
            parsing_struct& ps = p_node.ps[j];
            // 检查该节点是否需要建立分支
            if(!ps.essentially_empty()) { // 不可移动原子已插入 // FIXME ?!
                // 用main_branch创建residue并放入nr的flex向量中
                nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords)));
                postprocess_branch(nr, ps, c, nr.flex.back());
            }
        }
    }
    nr_update_matrixes(nr); // FIXME ?
    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

/**
 * @brief 从输入流解析配体PDBQT（dkoes的流版本）
 * @param in 输入流
 * @param nr 非刚性解析结果引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 支持从内存流解析配体，适用于动态生成的PDBQT数据
 */
void parse_pdbqt_ligand(std::istream& in, non_rigid_parsed& nr, context& c) {
    parsing_struct p;
    boost::optional<unsigned> torsdof;

    parse_pdbqt_aux(in, p, c, torsdof, false);

    if(p.atoms.empty()) 
        throw pdbqt_parse_error("No atoms in this ligand.");
    if(!torsdof)
        throw pdbqt_parse_error("Missing TORSDOF keyword.");

    postprocess_ligand(nr, p, c, unsigned(torsdof.get()));

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

/**
 * @brief 从文件解析配体PDBQT
 * @param name 文件路径
 * @param nr 非刚性解析结果引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 */
void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
    ifile in(name);
    parsing_struct p;
    boost::optional<unsigned> torsdof;

    parse_pdbqt_aux(in, p, c, torsdof, false);

    if(p.atoms.empty()) 
        throw pdbqt_parse_error("No atoms in this ligand.");
    if(!torsdof)
        throw pdbqt_parse_error("Missing TORSDOF keyword in this ligand.");

    postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // 奇怪的size_t -> unsigned编译器投诉

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

/**
 * @brief 从输入流解析柔性残基
 * @param in 输入流
 * @param p 解析结构引用
 * @param c 上下文引用
 */
void parse_pdbqt_residue(std::istream& in, parsing_struct& p, context& c) { 
    boost::optional<unsigned> dummy;
    parse_pdbqt_aux(in, p, c, dummy, true);
}

/**
 * @brief 解析柔性残基PDBQT文件
 * @param name 文件路径
 * @param nr 非刚性解析结果引用
 * @param c 上下文引用
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 柔性残基文件包含多个BEGIN_RES...END_RES块，每个块定义一个柔性残基
 */
void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
    ifile in(name);
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // 忽略空行
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "BEGIN_RES")) {
            parsing_struct p;
            parse_pdbqt_residue(in, p, c);
            postprocess_residue(nr, p, c);
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue PDBQT file. "
                                    "Use \"vina_split\" to split flex residues in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue.", str);
    }

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

/**
 * @brief 解析PDBQT分支
 * @param in 输入流
 * @param p 解析结构引用
 * @param c 上下文引用
 * @param from 起始原子编号
 * @param to 结束原子编号
 * @throws pdbqt_parse_error 解析失败时抛出异常
 * 
 * @note 递归解析BRANCH...ENDBRANCH之间的内容，支持嵌套分支结构
 */
void parse_pdbqt_branch(std::istream& in, parsing_struct& p, context& c, unsigned from, unsigned to) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // 忽略空行
        else if(starts_with(str, "WARNING")) {} // 忽略WARNING行 - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // 忽略REMARK行
        else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, str, p, c);
        else if(starts_with(str, "ENDBRANCH")) {
            unsigned first, second;
            parse_two_unsigneds(str, "ENDBRANCH", first, second);
            if(first != from || second != to) 
                throw pdbqt_parse_error("Inconsistent branch numbers.");
            if(!p.immobile_atom) 
                throw pdbqt_parse_error("Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch.");
            return;
        }
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
            parsed_atom a = parse_pdbqt_atom_string(str);
            if(a.number == to)
                p.immobile_atom = p.atoms.size();
            p.add(a, c);
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

//////////// new stuff //////////////////
/**
 * @struct pdbqt_initializer
 * @brief PDBQT解析初始化器，负责将解析后的数据转换为模型对象
 * 
 * 该结构体作为解析PDBQT文件后的数据处理中心，将原始解析数据转换为
 * AutoDock Vina内部使用的模型表示形式。
 */
struct pdbqt_initializer {

    atom_type::t atom_typing_used;  ///< 使用的原子类型系统
    model m;                        ///< 内部模型对象

    /**
     * @brief 构造函数
     * @param atype 原子类型系统
     */
    pdbqt_initializer(atom_type::t atype): atom_typing_used(atype), m(atype) {}

    /**
     * @brief 从刚性受体数据初始化模型
     * @param r 刚性受体结构数据
     * @note 将刚性受体的原子数据复制到模型的grid_atoms中，用于后续的网格计算
     */
    void initialize_from_rigid(const rigid& r) { // static really
        VINA_CHECK(m.grid_atoms.empty());
        m.grid_atoms = r.atoms;
    }

    /**
     * @brief 从非刚性解析数据初始化模型
     * @param nrp 非刚性解析结构数据
     * @param c 上下文信息
     * @param is_ligand 是否为配体（true为配体，false为柔性残基）
     * 
     * @note 该函数是模型初始化的核心，处理所有可移动和不可移动的原子：
     * - 复制配体和柔性残基信息
     * - 将可移动原子设置相对坐标，不可移动原子设置为零向量
     * - 初始化力计算所需的数据结构
     */
    void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
        VINA_CHECK(m.ligands.empty());
        VINA_CHECK(m.flex   .empty());

        // 复制配体和柔性残基数据
        m.ligands = nrp.ligands;
        m.flex    = nrp.flex;

        VINA_CHECK(m.atoms.empty());

        // 计算总原子数并预分配内存
        sz n = nrp.atoms.size() + nrp.inflex.size();
        m.atoms.reserve(n);
        m.coords.reserve(n);

        // 处理可移动原子：设置相对坐标用于旋转计算
        VINA_FOR_IN(i, nrp.atoms) {
            const movable_atom& a = nrp.atoms[i];
            atom b = static_cast<atom>(a);
            b.coords = a.relative_coords;  // 使用相对坐标
            m.atoms.push_back(b);
            m.coords.push_back(a.coords);  // 保存实际坐标
        }

        // 处理不可移动原子：坐标设为零向量避免混淆
        VINA_FOR_IN(i, nrp.inflex) {
            const atom& a = nrp.inflex[i];
            atom b = a;
            b.coords = zero_vec; // 避免混淆，这些坐标不会被访问
            m.atoms.push_back(b);
            m.coords.push_back(a.coords);
        }
        VINA_CHECK(m.coords.size() == n);

        // 初始化力计算相关数据结构
        m.minus_forces = m.coords;
        m.m_num_movable_atoms = nrp.atoms.size();

        // 根据类型设置上下文信息
        if(is_ligand) {
            VINA_CHECK(m.ligands.size() == 1);
            m.ligands.front().cont = c;  // 配体使用前端上下文
        }
        else
            m.flex_context = c;          // 柔性残基使用flex上下文
    }

    /**
     * @brief 使用移动性矩阵完成模型初始化
     * @param mobility 原子间移动性关系矩阵
     * @note 调用模型的initialize方法，设置原子间的键连关系和移动约束
     */
    void initialize(const distance_type_matrix& mobility) {
        m.initialize(mobility);
    }
};

/**
 * @brief 从文件解析配体PDBQT并返回模型
 * @param name PDBQT文件路径
 * @param atype 原子类型系统
 * @return model 解析后的配体模型
 * @throws parse_error 解析错误时抛出异常
 * 
 * @note 该函数是配体文件解析的主要入口点，完整流程包括：
 * 1. 解析PDBQT文件获取原始数据
 * 2. 初始化模型结构
 * 3. 设置原子移动性约束
 */
model parse_ligand_pdbqt_from_file(const std::string& name, atom_type::t atype) { // can throw parse_error
    non_rigid_parsed nrp;  // 非刚性解析结果
    context c;             // 上下文信息
    parse_pdbqt_ligand(make_path(name), nrp, c);  // 解析配体文件

    pdbqt_initializer tmp(atype);
    tmp.initialize_from_nrp(nrp, c, true);        // 作为配体初始化
    tmp.initialize(nrp.mobility_matrix());        // 设置移动性矩阵
    return tmp.m;
}

/**
 * @brief 从字符串解析配体PDBQT并返回模型
 * @param string_name 包含PDBQT内容的字符串
 * @param atype 原子类型系统
 * @return model 解析后的配体模型
 * @throws parse_error 解析错误时抛出异常
 */
model parse_ligand_pdbqt_from_string(const std::string& string_name, atom_type::t atype) { // can throw parse_error
    non_rigid_parsed nrp;
    context c;

    std::stringstream molstream(string_name);     // 将字符串转换为流
    parse_pdbqt_ligand(molstream, nrp, c);        // 从流解析配体

    pdbqt_initializer tmp(atype);
    tmp.initialize_from_nrp(nrp, c, true);
    tmp.initialize(nrp.mobility_matrix());
    return tmp.m;
}

/**
 * @brief 解析受体PDBQT文件（包含刚性部分和柔性残基）
 * @param rigid_name 刚性受体文件路径
 * @param flex_name 柔性残基文件路径
 * @param atype 原子类型系统
 * @return model 解析后的受体模型
 * 
 * @note 该函数处理受体的复杂结构：
 * - 刚性部分：固定不动的受体骨架，用于网格计算
 * - 柔性残基：可移动的侧链，参与对接优化
 * - 支持只有刚性部分或只有柔性部分的情况
 */
model parse_receptor_pdbqt(const std::string& rigid_name, const std::string& flex_name, atom_type::t atype) { 
    // 初始化数据结构
    rigid r;                    // 刚性受体数据
    non_rigid_parsed nrp;       // 非刚性解析数据
    context c;                  // 上下文信息
    pdbqt_initializer tmp(atype);

    // 解析刚性受体部分
    if (!rigid_name.empty()) {
        parse_pdbqt_rigid(make_path(rigid_name), r);
    }
    
    // 解析柔性残基部分
    if (!flex_name.empty()) {
        parse_pdbqt_flex(make_path(flex_name), nrp, c);
    }

    // 初始化刚性部分
    if (!rigid_name.empty()) {
        tmp.initialize_from_rigid(r);
        // 如果只有刚性部分，使用空的移动性矩阵
        if (flex_name.empty()) {
            distance_type_matrix mobility_matrix;
            tmp.initialize(mobility_matrix);
        }
    }

    // 初始化柔性部分
    if (!flex_name.empty()) {
        tmp.initialize_from_nrp(nrp, c, false);   // 不是配体，是柔性残基
        tmp.initialize(nrp.mobility_matrix());    // 使用柔性残基的移动性矩阵
    }

    return tmp.m;
}
