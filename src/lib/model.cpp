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

#include <random>

#include "model.h"
#include "file.h"
#include "curl.h"
#include "precalculate.h"

/**
 * @brief 深度优先计算原子树结构的原子索引范围
 * @tparam T 支持node和children成员的树结构类型
 * @param t 待计算的树结构节点
 * @return atom_range 包含所有子节点原子的范围
 */
template<typename T>
atom_range get_atom_range(const T& t) {
    atom_range tmp = t.node;
    VINA_FOR_IN(i, t.children) {
        atom_range r = get_atom_range(t.children[i]);  // 递归获取子节点范围
        // 扩展范围以包含子节点的所有原子
        if(tmp.begin > r.begin) tmp.begin = r.begin;
        if(tmp.end   < r.end  ) tmp.end   = r.end;
    }
    return tmp;  // 返回合并后的原子范围
}

/**
 * @brief 分子分支的几何度量结构
 */
struct branch_metrics {
    sz length;         ///< 分支长度：从根到最远叶子节点的距离
    sz corner2corner;  ///< 角到角距离：考虑分叉时两个最远端点间的距离
    
    /**
     * @brief 默认构造函数，初始化所有度量为0
     */
	branch_metrics() : length(0), corner2corner(0) {}
};

/**
 * @brief 递归计算分子分支的几何度量
 * 
 * 分析分子树结构，计算分支的长度和角到角距离，
 * 这些度量用于评估分子的空间分布特征
 * 
 * @tparam T 支持children成员的树结构类型
 * @param t 待分析的分子树节点
 * @return branch_metrics 包含长度和角到角距离的度量结果
 * 
 * @note 算法核心思想：
 *       1. length = 到最远子节点的最大距离
 *       2. corner2corner = 两个最长分支的长度之和，或单个最长分支长度
 * @note 对于有分叉的分子，corner2corner能更好地描述其空间延展
 */
template<typename T>
branch_metrics get_branch_metrics(const T& t) {
    branch_metrics tmp;  // 初始化结果结构
    
    if(!t.children.empty()) {  // 如果节点有子节点
        sz corner2corner_max = 0;  // 记录子树中的最大角到角距离
        szv lengths;  // 存储所有子分支的长度
        
        // 递归计算所有子节点的度量
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);

            // 更新子树中的最大角到角距离
			if(corner2corner_max < res.corner2corner)
				corner2corner_max = res.corner2corner;
            
            // 子分支长度加1（包含当前节点到子节点的连接）
			lengths.push_back(res.length + 1); // FIXME? weird compiler warning (sz -> unsigned)
		}
        
        // 按长度排序，便于找到最长的分支
		std::sort(lengths.begin(), lengths.end());

        // 分支长度是最长的子分支长度
		tmp.length = lengths.back();

        // 初始角到角距离等于最长分支长度
		tmp.corner2corner = tmp.length;
        // 如果有至少两个子分支，角到角距离是两个最长分支的和
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 2];  // 第二长的分支

        // 取当前计算值和子树最大值中的较大者
		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

/**
 * @brief 获取指定配体的最长分支长度
 * 
 * @param ligand_number 配体编号，必须小于配体总数
 * @return sz 配体最长分支的长度（原子数）
 * 
 * @note 最长分支长度反映了分子的线性延展程度
 */
sz model::ligand_longest_branch(sz ligand_number) const {
    return get_branch_metrics(ligands[ligand_number]).length;
}

/**
 * @brief 获取指定配体的角到角距离
 * 
 * @param ligand_number 配体编号，必须小于配体总数
 * @return sz 配体的角到角距离（原子数）
 * 
 * @note 角到角距离考虑了分子的分叉结构，能更好地描述分子的整体尺寸
 * @note 对于线性分子，角到角距离等于分支长度
 * @note 对于有分叉的分子，角到角距离通常大于单个分支长度
 */
sz model::ligand_length(sz ligand_number) const {
    return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

/**
 * @brief 调用get_atom_range设置配体的原子索引范围
 * @note 这个方法通常在配体初始化后调用
 */
void ligand::set_range() {
    atom_range tmp = get_atom_range(*this);  // 获取配体的原子范围
    begin = tmp.begin;  // 设置起始原子索引
    end   = tmp.end;    // 设置结束原子索引
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

/**
 * @brief 模型合并信息记录结构
 * 
 * 记录模型合并过程中需要的关键尺寸信息，用于索引转换和内存分配
 */
struct appender_info {
    sz grid_atoms_size;      ///< 网格原子数量
    sz m_num_movable_atoms;  ///< 可移动原子数量  
    sz atoms_size;           ///< 总原子数量

    /**
     * @brief 构造函数，从模型中提取尺寸信息
     * @param m 源模型对象
     */
	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

/**
 * @brief 模型合并执行器类
 * 
 * 负责将两个模型（a和b）合并为一个统一的模型。处理索引转换、坐标重组、
 * 相互作用对更新等复杂操作，确保合并后的模型保持物理和化学意义的正确性。
 * 
 * @note 合并策略：
 *       - 网格原子：a在前，b在后 (aaa...bbb...)
 *       - 可移动原子：交错排列 (aaa...bbb...aaa...)
 *       - 不可移动原子：a在前，b在后
 */
class appender {
	appender_info a_info;
	appender_info b_info;

    /**
     * @brief 计算网格原子的新索引
     * @param x 原始索引
     * @return 合并后的新索引
     * @note a的网格原子在b的网格原子之前拼接
     */
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
    bool is_a;  ///< 当前处理的是否为模型a的数据

    /**
     * @brief 构造函数，初始化合并器
     * @param a 第一个模型（目标模型）
     * @param b 第二个模型（源模型）
     */
	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

    /**
     * @brief 坐标索引转换操作符
     * @param x 原始坐标索引
     * @return 转换后的新索引
     * 
     * @note 转换规则：
     *       - 模型a：可移动原子不变，不可移动原子后移b的可移动原子数量
     *       - 模型b：可移动原子前移a的可移动原子数量，不可移动原子后移a的总原子数量
     */
	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
    /**
     * @brief 原子索引转换操作符
     * @param x 原始原子索引
     * @return 转换后的新原子索引
     * 
     * @note 区分网格原子和普通原子，分别应用不同的索引转换规则
     */
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
        if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);	// 网格原子使用网格索引转换
        else            tmp.i = operator()(tmp.i);		// 普通原子使用坐标索引转换
		return tmp;
	}

    // type-directed old -> new transformations
    /**
     * @brief 更新相互作用对的索引
     * @param ip 相互作用对引用
     */
	void update(interacting_pair& ip) const {
        ip.a = operator()(ip.a);  // 转换原子a的索引
        ip.b = operator()(ip.b);  // 转换原子b的索引
	}
    
    /**
     * @brief 更新向量（坐标和力）
     * @param v 向量引用
     * @note 向量本身不需要转换，只需转换其索引位置
     */
    void update(vec& v) const { // coordinates & forces - do nothing
    }
    
    /**
     * @brief 更新配体的所有相关索引
     * @param lig 配体引用
     * 
     * @note 更新内容包括：
     *       - 配体自身的原子范围
     *       - 内部相互作用对
     *       - 解析行信息
     */
    void update(ligand& lig) const {
        lig.transform(*this); // ligand as an atom_range subclass
        transform_ranges(lig, *this);  // 转换配体的原子范围
        VINA_FOR_IN(i, lig.pairs)
            this->update(lig.pairs[i]);  // 更新配体内部相互作用对
        VINA_FOR_IN(i, lig.cont)
            this->update(lig.cont[i]); // parsed_line update, below
    }
    
    /**
     * @brief 更新残基的原子范围
     * @param r 残基引用
     */
    void update(residue& r) const {
        transform_ranges(r, *this);
    }
    
    /**
     * @brief 更新解析行中的原子索引
     * @param p 解析行引用
     */
    void update(parsed_line& p) const {
        if(p.second)
            p.second = operator()(p.second.get());  // 转换可选的原子索引
    }
    
    /**
     * @brief 更新原子的键连接信息
     * @param a 原子引用
     */
    void update(atom& a) const {
        VINA_FOR_IN(i, a.bonds) {
            bond& b = a.bonds[i];
            b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
        }
    }

    /**
     * @brief 通用合并模板函数
     * @tparam T 容器元素类型
     * @param a 目标容器（第一个参数变成 aaaaaaaabbbbbbbbbbbbbbb）
     * @param b 源容器
     * 
     * @note 用于配体、柔性部分、柔性上下文、原子等的合并，也用于other_pairs
     */
    template<typename T>
    void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
        sz a_sz = a.size();  // 记录原始a的大小
        vector_append(a, b); // 直接拼接向量

        // 先更新原有的a部分数据
        is_a = true;
        VINA_FOR(i, a_sz)
            update(a[i]);

        // 再更新新添加的b部分数据
        is_a = false;
        VINA_RANGE(i, a_sz, a.size())
            update(a[i]);
    }

    /**
     * @brief 坐标特化的合并模板函数
     * @tparam T 坐标类型
     * @param a 目标坐标容器（第一个参数变成 aaaaaaaabbbbbbbbbaab）
     * @param b 源坐标容器
     * 
     * @note 坐标合并采用交错策略：
     *       1. a的可移动原子坐标保持在前
     *       2. b的可移动原子坐标插入到a的可移动原子后
     *       3. a的不可移动原子坐标
     *       4. b的不可移动原子坐标在最后
     */
    template<typename T>
    void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
        std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

        // 更新a的索引
        is_a = true;
        VINA_FOR_IN(i, a)
            update(a[i]);

        // 更新b副本的索引
        is_a = false;
        VINA_FOR_IN(i, b_copy)
            update(b_copy[i]);

        // 交错插入坐标
        typedef typename std::vector<T>::const_iterator const cci;
        cci b1 = b_copy.begin();                                    // b的开始
        cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;      // b的可移动原子结束
        cci b3 = b_copy.end();                                      // b的结束

        // 在a的可移动原子后插入b的可移动原子
        a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
        // 在最后插入b的不可移动原子
        a.insert(a.end()                                , b2 , b3);
    }
};

/**
 * @brief 将另一个模型合并到当前模型
 * @param m 要合并的源模型
 * 
 * @note 合并过程包括：
 *       1. 验证原子类型兼容性
 *       2. 合并各种相互作用对
 *       3. 重组坐标和力向量
 *       4. 合并配体、柔性部分等
 *       5. 添加新的分子间相互作用对
 *       6. 更新可移动原子计数
 */
void model::append(const model& m) {
    // 验证原子类型系统兼容性
    VINA_CHECK(atom_typing_used() == m.atom_typing_used());

    appender t(*this, m);  // 创建合并执行器

    // 合并各种相互作用对
    t.append(other_pairs, m.other_pairs);  // 分子内相互作用对
    t.append(inter_pairs, m.inter_pairs);  // 分子间相互作用对
    t.append(glue_pairs, m.glue_pairs);    // 胶合相互作用对

    // 验证坐标和力的一致性
    VINA_CHECK(minus_forces.size() == coords.size());
    VINA_CHECK(m.minus_forces.size() == m.coords.size());

    // 使用特化的坐标合并函数
    t.coords_append(coords, m.coords);
    t.coords_append(minus_forces, m.minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

    // 合并分子组件
    t.append(ligands, m.ligands);           // 配体
    t.append(flex, m.flex);                 // 柔性残基
    t.append(flex_context, m.flex_context); // 柔性上下文

    // 合并原子相关数据
    t.append(grid_atoms, m.grid_atoms);     // 网格原子
    t.coords_append(atoms, m.atoms);        // 原子（使用坐标合并策略）

    // 添加先前存在原子与新加入（很可能是配体）原子之间的相互作用对
    /* 相互作用类型：
    - flex     - ligand   : YES (inter_pairs)
    - flex_i   - flex_j   : YES (other_pairs) but append is used mostly for adding ligand
    - ligand_i - ligand_j : YES (inter_pairs)
    - macrocycle closure interactions: NO (1-2, 1-3, 1-4)
    */
    VINA_FOR(i, m_num_movable_atoms) {  // 遍历原有的可移动原子
        VINA_RANGE(j, m_num_movable_atoms, m_num_movable_atoms + m.m_num_movable_atoms) {  // 遍历新加入的可移动原子
            // 跳过大环闭合碰撞和不匹配的闭合虚原子
            if (is_closure_clash(i, j) || is_unmatched_closure_dummy(i, j)) continue;

            const atom& a = atoms[i];  // 原有原子
            const atom& b = atoms[j];  // 新原子

            // 获取原子类型
			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

            // 只处理有效的原子类型（排除氢原子等）
			if (t1 < n && t2 < n) {
                // 计算原子类型对索引
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
                
                // 根据原子类型和位置确定相互作用类型
				if (is_glue_pair(i, j)) {
                    // 胶合原子之间的相互作用
					glue_pairs.push_back(interacting_pair(type_pair_index, i, j)); // glue_i - glue_j
				} else if (is_atom_in_ligand(i) && is_atom_in_ligand(j)) {
                    // 配体内部相互作用
					inter_pairs.push_back(interacting_pair(type_pair_index, i, j)); // INTER: ligand_i - ligand_j
				} else if (is_atom_in_ligand(i) || is_atom_in_ligand(j)) {
                    // 配体-柔性部分相互作用
					inter_pairs.push_back(interacting_pair(type_pair_index, i, j)); // INTER: flex - ligand
                } else {
                    // 柔性部分内部相互作用
					other_pairs.push_back(interacting_pair(type_pair_index, i, j));
				}
			}
		}
	}

    // 更新可移动原子总数
	m_num_movable_atoms += m.m_num_movable_atoms;
}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////


atom_index model::sz_to_atom_index(sz i) const {
    if(i < grid_atoms.size()) return atom_index(i                    ,  true);	// 网格原子
    else                      return atom_index(i - grid_atoms.size(), false);	// 普通原子
}

/**
 * @brief 获取两个原子之间的距离类型
 * @param mobility 移动性矩阵
 * @param i 原子索引i
 * @param j 原子索引j
 * 
 * @note 首先判断i和j是否有刚体中的不可移动原子，若是则只能返回DISTANCE_FIXED或DISTANCE_VARIABLE
 * 若都不是刚体中的原子，则访问mobility矩阵获取它们之间的距离类型（在解析时确定）
 */
distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
    if(i.in_grid && j.in_grid) return DISTANCE_FIXED;  // 两个网格原子间距离固定
    if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;  // 网格原子与普通原子
    if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;  // 普通原子与网格原子
    
    // 两个普通原子的情况
    assert(!i.in_grid);
    assert(!j.in_grid);
    assert(i.i < atoms.size());
    assert(j.i < atoms.size());
    sz a = i.i;
    sz b = j.i;
    if(a == b) return DISTANCE_FIXED;  // 同一个原子
    return (a < b) ? mobility(a, b) : mobility(b, a);  // 查询移动性矩阵（确保索引顺序）
}


const vec& model::atom_coords(const atom_index& i) const {
    return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}


fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
    return vec_distance_sqr(atom_coords(a), atom_coords(b));
}


struct bond_less { // FIXME rm!?
	bool operator()(const bond& a, const bond& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};

bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const {
    fl r2 = distance_sqr_between(a, b);  // a和b间的距离平方
    
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
        if(a == c || b == c) continue;  // 跳过a和b本身
        
        // 检查c与a和b的距离类型
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
        
        // 如果c相对于a和b都不可移动，且距离都更近，则存在中间原子
		if(ac != DISTANCE_VARIABLE &&
			bc != DISTANCE_VARIABLE &&
			distance_sqr_between(a, c) < r2 &&
			distance_sqr_between(b, c) < r2)
			return true;
	}
	return false;
}

struct beads {
    fl radius_sqr;  ///< 珠子半径的平方
    std::vector<std::pair<vec, szv> > data;  ///< 珠子数据：<中心坐标, 原子索引列表>
    /**
     * @brief 构造函数
     * @param reserve_size 预留空间大小
     * @param radius_sqr_ 珠子半径的平方
     */
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }

    
	void add(sz index, const vec& coords) {
        // 寻找合适的现有珠子
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
                data[i].second.push_back(index);  // 加入现有珠子
                return;
            }
        }
        // 未找到合适珠子，创建新珠子
        std::pair<vec, szv> tmp;
        tmp.first = coords;
        tmp.second.push_back(index);
        data.push_back(tmp);
    }
};

/**
 * @brief 基于相对移动性、距离和共价长度分配化学键
 * @param mobility 移动性矩阵，定义原子间的相对移动关系
 */
void model::assign_bonds(const distance_type_matrix& mobility) {
    const fl bond_length_allowance_factor = 1.1;  // 给与键长相比默认值（两个原子共价半径之和）可能变化的范围
    sz n = grid_atoms.size() + atoms.size();  // 总原子数

    // 构建珠子数据结构并将所有原子添加到珠子中
    const fl bead_radius = 15;
    beads beads_instance(n, sqr(bead_radius));
    VINA_FOR(i, n) {
        atom_index i_atom_index = sz_to_atom_index(i);
        beads_instance.add(i, atom_coords(i_atom_index));
    }
    
    // 为每个原子分配键
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom& i_atom = get_atom(i_atom_index);

        // 设置共价半径与截断值
		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
        szv relevant_atoms;
        
        // 遍历所有珠子
		VINA_FOR_IN(b, beads_instance.data) {
            // 跳过距离过远的珠子
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
            
            // 遍历珠子内的原子
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom& j_atom = get_atom(j_atom_index);

                // 确定共价键长度和距离类型
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
                
                // 若距离可变则不太可能建立化学键
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
                    // 基于共价半径判断是否在键合范围内
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
        
        // 在相关原子中寻找真正的键合原子
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
            if(j <= i) continue; // 避免重复处理（每对原子只处理一次）
            
			atom_index j_atom_index = sz_to_atom_index(j);
			atom& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
            
            // 最终键分配判断：距离合适且无中间原子阻挡
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && 
               !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
                bool rotatable = (dt == DISTANCE_ROTOR);  // 是否为可旋转键
				fl length = std::sqrt(r2);
                
                // 为两个原子都添加键信息
				i_atom.bonds.push_back(bond(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond(i_atom_index, length, rotatable));
			}
		}
	}
}

/**
 * @brief 检查原子是否与HD（氢供体）原子键合
 * @param a 待检查的原子
 * @return 是否与HD原子键合
 */
bool model::bonded_to_HD(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}

/**
 * @brief 检查原子是否与杂原子（不是C/H）键合
 */
bool model::bonded_to_heteroatom(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

/**
 * @brief 根据元素类型、化学环境（键合情况）和AutoDock原子类型，
 * 为每个原子分配X-Score评分函数使用的原子类型。
 */
void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom& a = get_atom(ai);
        a.assign_el();  // 首先分配元素类型
        sz& x = a.xs;   // X-Score原子类型引用

        // 判断氢键性质
        bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score公式忽略SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

        // 根据元素类型分配X-Score原子类型
		switch(a.el) {
            case EL_TYPE_H    : break;  // 氢原子通常不参与评分
            case EL_TYPE_C    :{
                // 碳原子类型基于大环闭合标记和杂原子键合情况
                if     (a.ad == AD_TYPE_CG0){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG0 : XS_TYPE_C_H_CG0;}
                else if(a.ad == AD_TYPE_CG1){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG1 : XS_TYPE_C_H_CG1;}
                else if(a.ad == AD_TYPE_CG2){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG2 : XS_TYPE_C_H_CG2;}
                else if(a.ad == AD_TYPE_CG3){x = bonded_to_heteroatom(a) ? XS_TYPE_C_P_CG3 : XS_TYPE_C_H_CG3;}
                else                        {x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H;}
                break;
            }
            // 氮原子类型基于氢键供体/受体性质
            case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
            // 氧原子类型基于氢键供体/受体性质  
            case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
            case EL_TYPE_S    : x = XS_TYPE_S_P; break;   // 硫原子
            case EL_TYPE_P    : x = XS_TYPE_P_P; break;   // 磷原子
            case EL_TYPE_F    : x = XS_TYPE_F_H; break;   // 氟原子
            case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;  // 氯原子
            case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;  // 溴原子
            case EL_TYPE_I    : x = XS_TYPE_I_H; break;   // 碘原子
            case EL_TYPE_Si   : x = XS_TYPE_Si; break;    // 硅原子
            case EL_TYPE_At   : x = XS_TYPE_At; break;    // 砹原子
            case EL_TYPE_Met  : x = XS_TYPE_Met_D; break; // 金属原子
            case EL_TYPE_Dummy: {  // 虚原子（用于大环闭合）
                if      (a.ad == AD_TYPE_G0) x = XS_TYPE_G0;
                else if (a.ad == AD_TYPE_G1) x = XS_TYPE_G1;
                else if (a.ad == AD_TYPE_G2) x = XS_TYPE_G2;
                else if (a.ad == AD_TYPE_G3) x = XS_TYPE_G3;
                else if (a.ad == AD_TYPE_W)  x = XS_TYPE_SIZE; // W原子在XS类型中不存在
                else VINA_CHECK(false);
                break;
            }
            case EL_TYPE_SIZE : break;
            default: VINA_CHECK(false);  // 未知元素类型
        }
    }
}

/**
 * @brief 深度优先搜索与指定原子键合的所有原子（到指定深度）
 * 
 * @param a 起始原子索引
 * @param n 搜索深度（键的数量）
 * @param out 输出向量，存储找到的原子索引
 */
void model::bonded_to(sz a, sz n, szv& out) const {
    if(!has(out, a)) { // 避免重复添加
        out.push_back(a);
        if(n > 0)  // 还有搜索深度
            VINA_FOR_IN(i, atoms[a].bonds) {
                const bond& b = atoms[a].bonds[i];
                if(!b.connected_atom_index.in_grid)  // 只考虑非网格原子
                    bonded_to(b.connected_atom_index.i, n-1, out);  // 递归搜索
            }
    }
}

/**
 * @brief 查找与指定原子键合的所有原子
 * 
 * @param a 起始原子索引
 * @param n 搜索深度
 * @return 键合原子的索引向量
 */
szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}

/**
 * @brief 检查两个原子是否存在大环闭合碰撞 
 */
bool model::is_closure_clash(sz i, sz j) const {
	sz t1 = atoms[i].get(atom_type::AD);
	sz t2 = atoms[j].get(atom_type::AD);
    
    // 如果是G-CG配对，不算碰撞
    if ((t1==AD_TYPE_CG0 && t2==AD_TYPE_G0) || (t2==AD_TYPE_CG0 && t1==AD_TYPE_G0) ||
        (t1==AD_TYPE_CG1 && t2==AD_TYPE_G1) || (t2==AD_TYPE_CG1 && t1==AD_TYPE_G1) ||
        (t1==AD_TYPE_CG2 && t2==AD_TYPE_G2) || (t2==AD_TYPE_CG2 && t1==AD_TYPE_G2) ||
        (t1==AD_TYPE_CG3 && t2==AD_TYPE_G3) || (t2==AD_TYPE_CG3 && t1==AD_TYPE_G3))
			return false;
            
    // 查找两原子的直接邻居（1-2相互作用）
	szv neighbors_of_i = bonded_to(i, 1);
	szv neighbors_of_j = bonded_to(j, 1);
    
    // 检查i的邻居中是否有CG原子
    bool i_has_CG0 = false;
    bool i_has_CG1 = false;
    bool i_has_CG2 = false;
    bool i_has_CG3 = false;
    VINA_FOR_IN(index, neighbors_of_i) {
        sz type_i = atoms[neighbors_of_i[index]].get(atom_type::AD);
            if      (type_i==AD_TYPE_CG0) i_has_CG0=true;
            else if (type_i==AD_TYPE_CG1) i_has_CG1=true;
            else if (type_i==AD_TYPE_CG2) i_has_CG2=true;
            else if (type_i==AD_TYPE_CG3) i_has_CG3=true;
    }
    
    // 检查j的邻居中是否有与i邻居相同类型的CG原子
    VINA_FOR_IN(index, neighbors_of_j) {
        sz type_j = atoms[neighbors_of_j[index]].get(atom_type::AD);
        if ((type_j==AD_TYPE_CG0 && i_has_CG0) ||
            (type_j==AD_TYPE_CG1 && i_has_CG1) ||
            (type_j==AD_TYPE_CG2 && i_has_CG2) ||
            (type_j==AD_TYPE_CG3 && i_has_CG3))
            return true;  // 发现大环闭合碰撞
    }
    return false;
}

/**
 * @brief 检查两个原子是否为不匹配的闭合虚原子
 */
bool model::is_unmatched_closure_dummy(sz i, sz j) const {
    sz t1 = atoms[i].get(atom_type::AD);
    sz t2 = atoms[j].get(atom_type::AD);
    if ((t1==AD_TYPE_G0 && t2!=AD_TYPE_CG0) || (t2==AD_TYPE_G0 && t1!=AD_TYPE_CG0) ||
        (t1==AD_TYPE_G1 && t2!=AD_TYPE_CG1) || (t2==AD_TYPE_G1 && t1!=AD_TYPE_CG1) ||
        (t1==AD_TYPE_G2 && t2!=AD_TYPE_CG2) || (t2==AD_TYPE_G2 && t1!=AD_TYPE_CG2) ||
        (t1==AD_TYPE_G3 && t2!=AD_TYPE_CG3) || (t2==AD_TYPE_G3 && t1!=AD_TYPE_CG3))
        return true;
    else
        return false;
}

/**
 * @brief 检查是否为胶合原子对
 */
bool model::is_glue_pair(sz i, sz j) const {
    sz t1 = atoms[i].get(atom_type::AD);
    sz t2 = atoms[j].get(atom_type::AD);

    if ((t1 == AD_TYPE_CG0 && t2 == AD_TYPE_G0) || (t2 == AD_TYPE_CG0 && t1 == AD_TYPE_G0) ||
        (t1 == AD_TYPE_CG1 && t2 == AD_TYPE_G1) || (t2 == AD_TYPE_CG1 && t1 == AD_TYPE_G1) ||
        (t1 == AD_TYPE_CG2 && t2 == AD_TYPE_G2) || (t2 == AD_TYPE_CG2 && t1 == AD_TYPE_G2) ||
        (t1 == AD_TYPE_CG3 && t2 == AD_TYPE_G3) || (t2 == AD_TYPE_CG3 && t1 == AD_TYPE_G3))
        return true;
    else
        return false;
}

/**
 * @brief 初始化相互作用对
 * 
 * 建立分子系统中所有相关的非键相互作用对，这些相互作用对将用于能量计算。
 * 相互作用对的类型包括分子内、分子间和胶合相互作用。
 * 
 * @param mobility 移动性矩阵，定义原子间的相对移动关系
 * @note 算法流程：
 *       1. 为每个原子查找其1-4键合邻居
 *       2. 对每对原子，检查是否应建立相互作用对
 *       3. 排除键合邻居和特殊的大环结构
 *       4. 根据原子所属分子确定相互作用类型
 */
void model::initialize_pairs(const distance_type_matrix& mobility) {
    /* 相互作用类型：
    - ligand_i - ligand_i : 是（仅1-4）(ligand.pairs)
    - flex_i   - flex_i   : 是（仅1-4）(other_pairs)
    - flex_i   - flex_j   : 是 (other_pairs)
    - 大环内部闭合相互作用: 否（1-2, 1-3, 1-4）
    */
    VINA_FOR_IN(i, atoms) {
        sz i_lig = find_ligand(i);
        szv bonded_atoms = bonded_to(i, 3);     // BUG?排除了1-4及以下的相互作用

        VINA_RANGE(j, i + 1, atoms.size()) {
            if (mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
                if (is_closure_clash(i, j) || is_unmatched_closure_dummy(i, j)) continue;
                
                // 获取原子类型索引
                sz t1 = atoms[i].get(atom_typing_used());
                sz t2 = atoms[j].get(atom_typing_used());
                sz n  = num_atom_types(atom_typing_used());

                if (t1 < n && t2 < n) { 
                    // 两个原子类型在快速计算表中的索引
                    sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
                    interacting_pair ip(type_pair_index, i, j);
                    
                    if (is_glue_pair(i, j)) {
                        glue_pairs.push_back(ip);
                    } else if (i_lig < ligands.size() && find_ligand(j) == i_lig) {
                        ligands[i_lig].pairs.push_back(ip);
                    } else if (!is_atom_in_ligand(i) && !is_atom_in_ligand(j)) {
                        other_pairs.push_back(ip);
                    }
                }
            }
        }
    }
}

/**
 * @brief 初始化模型
 * @param mobility 移动性矩阵
 * 
 * @note 初始化顺序：
 *       1. 首先设置配体的原子索引范围（确定原子归属）
 *       2. 然后分配化学键（建立拓扑连接）
 *       3. 接着分配原子类型（基于键合环境）
 *       4. 最后初始化相互作用对（基于拓扑和类型）
 */
void model::initialize(const distance_type_matrix& mobility) {
    VINA_FOR_IN(i, ligands)
        ligands[i].set_range();  // 设置每个配体的原子范围
    assign_bonds(mobility);      // 分配化学键
    assign_types();             // 分配原子类型
    initialize_pairs(mobility); // 初始化相互作用对
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////


/**
 * @brief 计算所有配体的内部相互作用对总数
 * 
 * 遍历所有配体，累计每个配体内部的相互作用对数量
 * 
 * @return sz 内部相互作用对的总数
 * 
 * @note 内部相互作用对指配体分子内部原子间的非键相互作用（通常是1-4及以上）
 */
sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();  // 累加每个配体的相互作用对数量
	return tmp;
}

/**
 * @brief 获取可移动原子的所有原子类型
 * @param atom_typing_used_ 使用的原子类型系统
 * @return szv 包含所有可移动原子类型的向量
 */
szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
    sz n = num_atom_types(atom_typing_used_);  // 获取指定类型系统的原子类型总数
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a = atoms[i];
        sz t = a.get(atom_typing_used_);  // 获取原子的类型索引
        if(t < n && !has(tmp, t))  // 类型有效且未添加过
			tmp.push_back(t);
	}
	return tmp;
}

/**
 * @brief 获取模型的构象大小信息
 * 
 * 统计配体和柔性部分的可旋转键（扭转角）数量
 * 
 * @return conf_size 包含配体和柔性部分自由度的结构
 * 
 * @note 构象大小决定了优化过程中需要搜索的维度数量
 */
conf_size model::get_size() const {
    conf_size tmp;
    tmp.ligands = ligands.count_torsions();  // 统计所有配体的扭转角总数
    tmp.flex    = flex   .count_torsions();  // 统计所有柔性残基的扭转角总数
    return tmp;
}

/**
 * @brief 获取初始构象
 * 
 * 创建一个扭转角为0、取向为单位矩阵、配体位置为当前位置的初始构象
 * 
 * @return conf 初始构象参数
 * 
 * @note 用于优化算法的起始点，确保从合理的初始状态开始搜索
 */
conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();  // 将所有扭转角设为0，取向设为单位矩阵
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();  // 设置配体当前位置
	return tmp;
}

/**
 * @brief 获取配体坐标（向量版本）
 * 
 * 返回第一个配体的所有原子坐标
 * 
 * @return vecv 配体原子坐标向量
 * 
 * @note FIXME 标记表明此函数可能需要重构或移除
 * @note 当前实现限制只能处理单个配体的情况
 */
vecv model::get_ligand_coords() const { // FIXME rm
    VINA_CHECK(ligands.size() == 1);  // 验证只有一个配体
	vecv tmp;
	const ligand &lig = ligands.front();
    VINA_RANGE(i, lig.begin, lig.end)  // 遍历配体的原子范围
        tmp.push_back(coords[i]);
    return tmp;
}

/**
 * @brief 获取配体坐标（双精度版本）
 * 
 * 将配体坐标转换为C++标准库格式，用于与外部C++代码交互
 * 
 * @return std::vector<double> 包含所有配体原子坐标的一维数组
 * 
 * @note 坐标按x,y,z顺序连续存储：[x1,y1,z1,x2,y2,z2,...]
 * @note 用于将坐标数据从C++内部格式导出到外部世界
 */
std::vector<double> model::get_ligand_coords() {
    // Way to get coordinates out of the C++ world
    VINA_CHECK(ligands.size() == 1);  // 验证只有一个配体
    std::vector<double> tmp;
    const ligand &lig = ligands.front();
    VINA_RANGE(i, lig.begin, lig.end) {
        tmp.push_back(coords[i][0]);  // x坐标
        tmp.push_back(coords[i][1]);  // y坐标
        tmp.push_back(coords[i][2]);  // z坐标
    }
    return tmp;
}

/**
 * @brief 获取重原子的可移动坐标
 * 
 * 提取所有可移动的重原子（非氢原子）坐标
 * 
 * @return vecv 重原子坐标向量
 * 
 * @note FIXME 标记表明此函数可能需要移动到其他位置
 * @note 重原子在结构分析和能量计算中更为重要
 */
vecv model::get_heavy_atom_movable_coords() const { // FIXME mv
    vecv tmp;

    VINA_FOR(i, num_movable_atoms()) {
        if (atoms[i].el != EL_TYPE_H)  // 排除氢原子
            tmp.push_back(coords[i]);
    }
    return tmp;
}

/**
 * @brief 通过原子索引查找原子所属的配体编号，如果不属于任何配体则返回配体总数
 */
sz model::find_ligand(sz a) const {
    VINA_FOR_IN(i, ligands) {
        if(a >= ligands[i].begin && a < ligands[i].end)
            return i;
    }
    return ligands.size();  // 返回无效索引表示不属于任何配体
}

/**
 * @brief 判断原子是否属于配体
 * 
 * 检查给定原子是否属于任何一个配体
 * 
 * @param a 原子索引
 * @return bool 是否属于配体
 */
bool model::is_atom_in_ligand(sz a) const {
    VINA_FOR_IN(i, ligands) {
        if (a >= ligands[i].begin && a < ligands[i].end)
            return true;
    }
    return false;
}

/**
 * @brief 判断原子是否可移动
 * 
 * 基于原子索引判断其是否为可移动原子
 * 
 * @param a 原子索引
 * @return bool 是否可移动
 * 
 * @note 可移动原子索引范围：[0, m_num_movable_atoms)
 */
bool model::is_movable_atom(sz a) const {
    if (a < num_movable_atoms())
        return true;
    else
        return false;
}

/**
 * @brief 计算可移动原子的几何中心
 * 
 * 计算所有可移动原子的质心坐标
 * 
 * @return std::vector<double> 几何中心坐标[x, y, z]
 * 
 * @note 几何中心 = Σ(坐标) / 原子数，每个原子权重相等
 */
std::vector<double> model::center() const {
    std::vector<double> center(3, 0);  // 初始化为[0, 0, 0]

    VINA_FOR(i, num_movable_atoms()) {
        center[0] += coords[i][0];  // 累加x坐标
        center[1] += coords[i][1];  // 累加y坐标
        center[2] += coords[i][2];  // 累加z坐标
    }

    // 计算平均值得到几何中心
    center[0] /= num_movable_atoms();
    center[1] /= num_movable_atoms();
    center[2] /= num_movable_atoms();

    return center;
}

/**
 * @brief 将坐标写入字符串的指定位置
 * 
 * 在PDBQT格式字符串的指定位置写入格式化的坐标值
 * 
 * @param i 起始位置（1基索引）
 * @param x 要写入的坐标值
 * @param str 目标字符串（会被修改）
 * 
 * @note PDBQT格式要求坐标占用8个字符，精度为3位小数
 * @note 使用1基索引，内部转换为0基索引
 */
void string_write_coord(sz i, fl x, std::string& str) {
    VINA_CHECK(i > 0);
    --i;  // 转换为0基索引
    std::ostringstream out;
    out.setf(std::ios::fixed, std::ios::floatfield);  // 固定小数点格式
    out.setf(std::ios::showpoint);  // 显示小数点
    out << std::setw(8) << std::setprecision(3) << x;  // 8字符宽度，3位精度
    VINA_CHECK(out.str().size() == 8);  // 验证格式化结果长度
    VINA_CHECK(str.size() > i + 8);     // 验证目标字符串有足够空间
    VINA_FOR(j, 8)
        str[i+j] = out.str()[j];  // 逐字符复制到目标位置
}

/**
 * @brief 将坐标转换为PDBQT格式字符串
 * 
 * 在现有PDBQT行字符串中更新坐标信息
 * 
 * @param coords 新的坐标值
 * @param str 原始PDBQT行字符串
 * @return std::string 更新坐标后的PDBQT字符串
 * 
 * @note PDBQT格式中坐标位置：
 *       - X坐标：第31-38列
 *       - Y坐标：第39-46列  
 *       - Z坐标：第47-54列
 */
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) {
    std::string tmp(str);
    string_write_coord(31, coords[0], tmp);  // 写入X坐标
    string_write_coord(39, coords[1], tmp);  // 写入Y坐标
    string_write_coord(47, coords[2], tmp);  // 写入Z坐标
    return tmp;
}

/**
 * @brief 将上下文写入文件输出流
 * 
 * 将PDBQT上下文（包含原始行和可选的原子索引）写入文件
 * 
 * @param c 要写入的上下文
 * @param out 输出文件流
 * 
 * @note 对于包含原子索引的行，会用当前坐标更新其坐标信息
 * @note 对于普通行，直接写入原始字符串
 */
void model::write_context(const context& c, ofile& out) const {
    verify_bond_lengths();  // 验证键长的一致性
    VINA_FOR_IN(i, c) {
        const std::string& str = c[i].first;
        if(c[i].second) {
            // 包含原子索引的行：更新坐标信息
            out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
        }
        else
            out << str << '\n';  // 普通行：直接输出
    }
}

/**
 * @brief 将上下文写入字符串流
 * 
 * 将PDBQT上下文写入字符串流，用于生成字符串格式的输出
 * 
 * @param c 要写入的上下文
 * @param out 输出字符串流
 */
void model::write_context(const context &c, std::ostringstream& out) const {
    verify_bond_lengths();  // 验证键长的一致性

    VINA_FOR_IN(i, c) {
        const std::string &str = c[i].first;
        if (c[i].second)
            // 包含原子索引的行：更新坐标信息
            out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
        else
            out << str << '\n';  // 普通行：直接输出
    }
}

/**
 * @brief 写入模型到字符串
 * 
 * 生成包含模型信息的PDBQT格式字符串
 * 
 * @param model_number 模型编号
 * @param remark 备注信息
 * @return std::string 完整的模型字符串
 * 
 * @note 输出格式：
 *       MODEL <编号>
 *       <备注信息>
 *       <配体信息>
 *       <柔性部分信息>（如果存在）
 *       ENDMDL
 */
std::string model::write_model(sz model_number, const std::string &remark) {
	std::ostringstream out;

	out << "MODEL " << model_number << '\n';
	out << remark;

    // 写入所有配体的上下文
	VINA_FOR_IN(i, ligands)
		write_context(ligands[i].cont, out);
    
    // 如果有柔性部分，写入其上下文
	if (num_flex() > 0) // otherwise remark is written in vain
		write_context(flex_context, out);

	out << "ENDMDL\n";

	return out.str();
}

/**
 * @brief 设置模型构象
 * 
 * 根据给定的构象参数更新原子坐标
 * 
 * @param c 构象参数，包含配体和柔性部分的姿态信息
 * 
 * @note 该方法会同时更新配体和柔性部分的原子坐标
 */
void model::set         (const conf& c) {
    ligands.set_conf(atoms, coords, c.ligands);  // 设置配体构象
    flex   .set_conf(atoms, coords, c.flex);     // 设置柔性部分构象
}

/**
 * @brief 计算指定配体的回转半径
 * 
 * 回转半径是描述分子空间分布的重要几何参数，定义为所有重原子到分子质心距离的均方根
 * 
 * @param ligand_number 配体编号，必须小于配体总数
 * @return 回转半径值(Å)，如果没有重原子则返回0
 * 
 * @note 只考虑重原子(非氢原子)，使用配体节点的几何中心作为参考点
 * @note 计算公式: Rg = sqrt(Σ(ri - r_center)² / N)，其中ri是原子坐标，r_center是质心
 */
fl model::gyration_radius(sz ligand_number) const {
    VINA_CHECK(ligand_number < ligands.size());  // 验证配体编号有效性
    const ligand& lig = ligands[ligand_number];  // 获取指定配体
    fl acc = 0;         // 累积距离平方和
    unsigned counter = 0;  // 重原子计数器
    
    // 遍历配体中的所有原子
    VINA_RANGE(i, lig.begin, lig.end) {
        if(atoms[i].el != EL_TYPE_H) { // 只处理重原子(非氢原子)
            // 计算原子到配体几何中心的距离平方并累加
            acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
            ++counter;  // 增加重原子计数
        }
    }
    // 返回均方根距离，如果没有重原子则返回0
    return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


/**
 * @brief 评估原子间相互作用对的总能量
 * 
 * 遍历所有相互作用对，计算在截断距离内的原子对的相互作用能量
 * 
 * @param p 预计算的原子间相互作用参数
 * @param v 能量调节参数，用于curl函数
 * @param pairs 需要评估的相互作用对列表
 * @param coords 所有原子的坐标数组
 * @param with_max_cutoff 是否使用最大截断距离(默认false)
 * @return 总的相互作用能量
 * 
 * @note 使用截断距离优化：只计算距离小于cutoff的原子对
 * @note curl函数用于能量的平滑截断，避免能量突变
 */
fl eval_interacting_pairs(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, const bool with_max_cutoff) { // clean up
    fl e = 0;  // 总能量初始化
    fl cutoff_sqr = p.cutoff_sqr();  // 获取截断距离的平方

    // 根据参数选择使用普通截断距离还是最大截断距离
	if (with_max_cutoff) {
		cutoff_sqr = p.max_cutoff_sqr();
	}

    // 遍历所有相互作用对
	VINA_FOR_IN(i, pairs) {
        const interacting_pair& ip = pairs[i];  // 获取当前相互作用对
        // 计算两原子间的距离平方
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);

        if(r2 < cutoff_sqr) {  // 只处理在截断距离内的原子对
            // 使用预计算参数快速评估相互作用能量
			fl tmp = p.eval_fast(ip.a, ip.b, r2);
            curl(tmp, v);  // 应用能量平滑截断
            e += tmp;      // 累加到总能量
		}
	}
	return e;
}

/**
 * @brief 评估原子间相互作用对的能量和导数(力)
 * 
 * 不仅计算相互作用能量，还计算作用在每个原子上的力，用于分子动力学优化
 * 
 * @param p 预计算的原子间相互作用参数
 * @param v 能量调节参数，用于curl函数
 * @param pairs 需要评估的相互作用对列表
 * @param coords 所有原子的坐标数组
 * @param forces 输出参数：各原子受到的力向量(会累加到现有值)
 * @param with_max_cutoff 是否使用最大截断距离(默认false)
 * @return 总的相互作用能量
 * 
 * @note 力的计算基于能量对坐标的负梯度: F = -∇E
 * @note 使用牛顿第三定律：作用力和反作用力大小相等方向相反
 * @note curl函数同时处理能量和力的平滑截断
 */
fl eval_interacting_pairs_deriv(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces, const bool with_max_cutoff) { // adds to forces  // clean up
    fl e = 0;  // 总能量初始化
    fl cutoff_sqr = p.cutoff_sqr();  // 获取截断距离的平方

    // 根据参数选择使用普通截断距离还是最大截断距离
    if (with_max_cutoff) {
		cutoff_sqr = p.max_cutoff_sqr();
	}

    // 遍历所有相互作用对
	VINA_FOR_IN(i, pairs) {
        const interacting_pair& ip = pairs[i];  // 获取当前相互作用对
        vec r = coords[ip.b] - coords[ip.a];    // 计算从原子a指向原子b的向量
        fl r2 = sqr(r);  // 计算距离向量的模长平方
        
        if(r2 < cutoff_sqr) {  // 只处理在截断距离内的原子对
            // 计算能量和其对距离平方的导数
            pr tmp = p.eval_deriv(ip.a, ip.b, r2);  // tmp.first是能量，tmp.second是dE/dr²
            // 计算作用力向量：F = -dE/dr = -(dE/dr²) * (dr²/dr) = -(dE/dr²) * 2r
			vec force;
            force = tmp.second * r;  // 力向量 = 导数 × 位置向量
            
            curl(tmp.first, force, v);  // 对能量和力应用平滑截断
            e += tmp.first;  // 累加能量

            // 根据牛顿第三定律分配力：
            // FIXME 如果使用硬截断，这种方式效率较低
            forces[ip.a] -= force;  // 原子a受到指向原子b的力(负号因为力的方向)
            forces[ip.b] += force;  // 原子b受到来自原子a的反作用力
                                   // 注意：我们可以忽略不可移动原子上的力
		}
	}
	return e;
}
/**
 * @brief 评估外部相互作用对能量（分子内其他相互作用）
 * 
 * 计算柔性残基间、柔性残基内部等其他类型的相互作用能量
 * 
 * @param p 预计算的原子间相互作用参数
 * @param v 能量调节参数向量，v[2]用于other_pairs的能量调节
 * @return 其他相互作用对的总能量
 * 
 * @note 主要用于计算受体柔性部分的内部相互作用
 */
fl model::evalo(const precalculate_byatom& p, const vec& v) const { // clean up
    fl e = eval_interacting_pairs(p, v[2], other_pairs, coords);
    return e;
}

/**
 * @brief 评估分子间相互作用能量
 * 
 * 计算配体与柔性残基间、不同配体间的相互作用能量
 * 
 * @param p 预计算的原子间相互作用参数
 * @param v 能量调节参数向量，v[2]用于inter_pairs的能量调节
 * @return 分子间相互作用的总能量
 * 
 * @note 这是配体与受体间主要相互作用能量的重要组成部分
 */
fl model::eval_inter(const precalculate_byatom& p, const vec& v) const { // clean up
    fl e = eval_interacting_pairs(p, v[2], inter_pairs, coords);
    return e;
}

/**
 * @brief 评估配体内部相互作用能量
 * 
 * 计算所有配体分子内部的非键相互作用能量（通常是1-4及以上相互作用）
 * 
 * @param p 预计算的原子间相互作用参数  
 * @param v 能量调节参数向量，v[0]用于配体内部相互作用的能量调节
 * @return 所有配体内部相互作用的总能量
 * 
 * @note 遍历每个配体，累加其内部相互作用对的能量贡献
 */
fl model::evali(const precalculate_byatom& p, const vec& v) const { // clean up
    fl e = 0;
    VINA_FOR_IN(i, ligands) 
        e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // probably might was well use coords here
    return e;
}

/**
 * @brief 计算总能量及其对坐标的导数（用于优化）
 * 
 * 计算系统总能量的同时，计算能量对各原子坐标的梯度，
 * 为优化算法提供力的信息
 * 
 * @param p 预计算的原子间相互作用参数
 * @param ig 交互网格，用于配体-受体刚性部分的相互作用
 * @param v 能量调节参数向量
 * @param g 输出参数：梯度变化，包含配体和柔性部分的力信息
 * @return 系统总能量
 * 
 * @note 能量组成包括：
 *       - 配体与网格的相互作用
 *       - 配体内部相互作用  
 *       - 配体间相互作用
 *       - 柔性部分相互作用
 *       - 胶合原子相互作用
 */
fl model::eval_deriv(const precalculate_byatom& p, const igrid& ig, const vec& v, change& g) { // clean up
    // INTER ligand - grid (配体与网格的相互作用)
    fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex

    // INTRA ligand_i - ligand_i (配体内部相互作用)
    VINA_FOR_IN(i, ligands)
        e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces

    // INTER ligand_i - ligand_j and ligand_i - flex_i (配体间及配体-柔性相互作用)
    if (!inter_pairs.empty()) 
        e += eval_interacting_pairs_deriv(p, v[2], inter_pairs, coords, minus_forces); // adds to minus_forces
    // INTRA flex_i - flex_i and flex_i - flex_j (柔性部分内部及柔性部分间相互作用)
    if (!other_pairs.empty())
        e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
    // glue_i - glue_i and glue_i - glue_j (胶合原子相互作用，用于大环闭合)
    if (!glue_pairs.empty())
        e += eval_interacting_pairs_deriv(p, v[2], glue_pairs, coords, minus_forces, true); // adds to minus_forces

    // calculate derivatives (计算最终的导数)
    ligands.derivative(coords, minus_forces, g.ligands);  // 计算配体的导数
    flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored，计算柔性部分的导数
    return e;
}

/**
 * @brief 计算分子内相互作用能量
 * 
 * 计算系统的分子内能量，包括配体内部、柔性部分内部以及柔性部分与刚性受体的相互作用
 * 
 * @param p 预计算的原子间相互作用参数
 * @param ig 交互网格，用于柔性-刚性相互作用
 * @param v 能量调节参数向量
 * @return 分子内总能量
 * 
 * @note 与总能量计算不同，这里只计算分子内部的能量贡献，
 *       不包括配体与刚性受体网格的相互作用
 */
fl model::eval_intramolecular(const precalculate_byatom& p, const igrid& ig, const vec& v) {
    fl e = 0;
    const fl cutoff_sqr = p.cutoff_sqr();

    // internal for each ligand (每个配体的内部相互作用)
    VINA_FOR_IN(i, ligands)
        e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

    // flex - rigid (柔性部分与刚性部分的相互作用)
    e += ig.eval_intra(*this, v[1]);

    // flex_i - flex_i and flex_i - flex_j (柔性部分的内部和相互间的相互作用)
    VINA_FOR_IN(i, other_pairs) {
        const interacting_pair& pair = other_pairs[i];
        fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
        if (r2 < cutoff_sqr) {
            fl this_e = p.eval_fast(pair.a, pair.b, r2);
            curl(this_e, v[2]);  // 应用能量平滑截断
            e += this_e;
        }
    }

    return e;
}

/**
 * @brief 计算两个模型间的非对称RMSD下界
 * 
 * 通过寻找x中每个重原子在y中最接近的同类型原子来计算RMSD下界
 * 
 * @param x 第一个模型
 * @param y 第二个模型
 * @return 非对称RMSD下界值（Å）
 * 
 * @note 算法为每个x中的重原子找到y中距离最近的同类型重原子，
 *       然后计算这些最小距离的均方根
 * @note 实际上是静态方法，但为了访问私有成员而声明为成员函数
 */
fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
    sz n = x.m_num_movable_atoms; 
    VINA_CHECK(n == y.m_num_movable_atoms);
    fl sum = 0;
    unsigned counter = 0;
    VINA_FOR(i, n) {
        const atom& a =   x.atoms[i];
        if(a.el != EL_TYPE_H) {  // 只考虑重原子（非氢原子）
            fl r2 = max_fl;
            VINA_FOR(j, n) {
                const atom& b = y.atoms[j];
                if(a.same_element(b) && !b.is_hydrogen()) {  // 寻找同类型的重原子
                    fl this_r2 = vec_distance_sqr(x.coords[i], 
                                                  y.coords[j]);
                    if(this_r2 < r2)
                        r2 = this_r2;  // 更新最小距离平方
                }
            }
            assert(not_max(r2));
            sum += r2;      // 累加距离平方
            ++counter;      // 计数重原子数量
        }
    }
    return (counter == 0) ? 0 : std::sqrt(sum / counter);  // 返回均方根距离
}

/**
 * @brief 计算两个模型间的RMSD下界
 * 
 * 取两个方向非对称RMSD下界的最大值作为最终下界
 * 
 * @param m 比较的模型
 * @return RMSD下界值（Å）
 * 
 * @note 由于非对称计算的方向性，取两个方向的最大值可获得更紧的下界
 */
fl model::rmsd_lower_bound(const model& m) const {
    return (std::max)(rmsd_lower_bound_asymmetric(*this, m), rmsd_lower_bound_asymmetric(m, *this));
}

/**
 * @brief 计算两个模型间的RMSD上界
 * 
 * 通过对应原子间的直接距离计算RMSD上界
 * 
 * @param m 比较的模型
 * @return RMSD上界值（Å）
 * 
 * @note 假设两个模型的原子顺序相同，直接计算对应原子间的距离
 * @note 只考虑重原子的贡献
 */
fl model::rmsd_upper_bound(const model& m) const {
    VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
    fl sum = 0;
    unsigned counter = 0;
    VINA_FOR(i, m_num_movable_atoms) {
        const atom& a =   atoms[i];
        const atom& b = m.atoms[i];
        assert(a.ad == b.ad);  // 验证原子类型相同
        assert(a.xs == b.xs);
        if(a.el != EL_TYPE_H) {  // 只考虑重原子
            sum += vec_distance_sqr(coords[i], m.coords[i]);  // 直接计算对应原子间距离平方
            ++counter;
        }
    }
    return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/**
 * @brief 计算配体部分的RMSD上界
 * 
 * 只计算配体原子的RMSD，忽略柔性受体部分
 * 
 * @param m 比较的模型
 * @return 配体RMSD上界值（Å）
 * 
 * @note 用于专门评估配体构象的变化，不受受体柔性部分影响
 */
fl model::rmsd_ligands_upper_bound(const model& m) const {
    VINA_CHECK(ligands.size() == m.ligands.size());
    fl sum = 0;
    unsigned counter = 0;
    VINA_FOR_IN(ligand_i, ligands) {
        const ligand&   lig =   ligands[ligand_i];
        const ligand& m_lig = m.ligands[ligand_i];
        VINA_CHECK(lig.begin == m_lig.begin);
        VINA_CHECK(lig.end   == m_lig.end);
        VINA_RANGE(i, lig.begin, lig.end) {  // 遍历配体的原子范围
            const atom& a =   atoms[i];
            const atom& b = m.atoms[i];
            assert(a.ad == b.ad);
            assert(a.xs == b.xs);
            if(a.el != EL_TYPE_H) {  // 只考虑重原子
                sum += vec_distance_sqr(coords[i], m.coords[i]);
                ++counter;
            }
        }
    }
    return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/**
 * @brief 验证模型中所有键的长度
 * 
 * 检查每个原子的每个键，验证计算的距离与存储的键长是否一致
 * 
 * @note 用于调试和验证模型的结构完整性
 * @note 如果发现不一致会显示详细信息并断言失败
 */
void model::verify_bond_lengths() const {
    VINA_FOR(i, grid_atoms.size() + atoms.size()) {
        const atom_index ai = sz_to_atom_index(i);
        const atom& a = get_atom(ai);
        VINA_FOR_IN(j, a.bonds) {
            const bond& b = a.bonds[j];
            fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
            bool ok = eq(d, b.length);  // 检查计算距离与存储距离是否相等
            if(!ok) {
                VINA_SHOW(d);           // 显示计算的距离
                VINA_SHOW(b.length);    // 显示存储的距离
            }
            VINA_CHECK(ok);  // 断言检查
        }
    }
}

/**
 * @brief 检查配体内部相互作用对的有效性
 * 
 * 验证每个配体内部的相互作用对中的原子索引是否在配体的原子范围内
 * 
 * @note 用于调试，确保相互作用对的索引正确性
 */
void model::check_ligand_internal_pairs() const {
    VINA_FOR_IN(i, ligands) {
        const ligand& lig = ligands[i];
        const interacting_pairs& pairs = lig.pairs;
        VINA_FOR_IN(j, pairs) {
            const interacting_pair& ip = pairs[j];
            VINA_CHECK(ip.a >= lig.begin);  // 验证原子a在配体范围内
            VINA_CHECK(ip.b  < lig.end);    // 验证原子b在配体范围内
        }
    }
}

/**
 * @brief 显示模型的基本信息
 * 
 * 输出模型的统计信息，包括原子类型、原子数量、相互作用对数量等
 * 
 * @note 用于调试和模型分析
 */
void model::about() const {
    VINA_SHOW(atom_typing_used());    // 显示使用的原子类型系统
    VINA_SHOW(num_movable_atoms());   // 显示可移动原子数量
    VINA_SHOW(num_internal_pairs());  // 显示内部相互作用对数量
    VINA_SHOW(num_other_pairs());     // 显示其他相互作用对数量
    VINA_SHOW(num_ligands());         // 显示配体数量
    VINA_SHOW(num_flex());            // 显示柔性残基数量
}

/**
 * @brief 显示模型中所有相互作用对的详细信息
 * 
 * 输出各类相互作用对的原子索引、类型和坐标信息，
 * 包括分子间、分子内配体、分子内柔性和胶合相互作用对
 * 
 * @note 用于调试相互作用对的生成和分类是否正确
 */
void model::show_pairs() const {

    std::cout << "INTER PAIRS\n";  // 分子间相互作用对
    interacting_pairs inter_pairs = get_inter_pairs();
    VINA_FOR_IN(i, inter_pairs) {
        const interacting_pair &ip = inter_pairs[i];

        if (is_atom_in_ligand(ip.a)) {
            sz lig_i = find_ligand(ip.a);
            std::cout << "LIGAND (" << lig_i << ") : ";
        } else {
            std::cout << "  FLEX     : ";
        }

        if (is_atom_in_ligand(ip.b)) {
            sz lig_i = find_ligand(ip.b);
            std::cout << "LIGAND (" << lig_i << ") ";
        } else {
            std::cout << "  FLEX     ";
        }

        std::cout << " - " << ip.a << " : " << ip.b
                  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
                  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
                  << "\n";
    }

    std::cout << "INTRA LIG PAIRS\n";  // 配体内部相互作用对
    VINA_FOR(i, num_ligands()) {
        ligand lig = get_ligand(i);
        VINA_FOR_IN(j, lig.pairs) {
            const interacting_pair &ip = lig.pairs[j];
            sz lig_i = find_ligand(ip.a);
            std::cout << "LIGAND (" << lig_i << ") ";
            std::cout << " - " << ip.a << " : " << ip.b
                      << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
                      << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
                      << "\n";
        }
    }

    std::cout << "INTRA FLEX PAIRS\n";  // 柔性部分内部相互作用对
    interacting_pairs other_pairs = get_other_pairs();
    VINA_FOR_IN(i, other_pairs) {
        const interacting_pair& ip = other_pairs[i];
        std::cout << "FLEX       ";
        std::cout << " - " << ip.a << " : " << ip.b
                  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
                  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
                  << "\n";
    }

    std::cout << "GLUE - GLUE PAIRS\n";  // 胶合相互作用对（大环闭合）
    interacting_pairs glue_pairs = get_glue_pairs();
    VINA_FOR_IN(i, glue_pairs) {
        const interacting_pair& ip = glue_pairs[i];
        std::cout << "FLEX       ";
        std::cout << " - " << ip.a << " : " << ip.b
                  << " - " << get_coords(ip.a)[0] << " " << get_coords(ip.a)[1] << " " << get_coords(ip.a)[2]
                  << " - " << get_coords(ip.b)[0] << " " << get_coords(ip.b)[1] << " " << get_coords(ip.b)[2]
                  << "\n";
    }
}

/**
 * @brief 显示模型中所有原子的详细信息
 * 
 * 输出每个原子的移动性、坐标、原子类型等信息
 * 
 * @note 用于调试原子分配和坐标检查
 */
void model::show_atoms() const {
    std::cout << "ATOM INFORMATION\n";
    VINA_FOR_IN(i, atoms) {
        const atom &a = atoms[i];
        if (i < num_movable_atoms()) {
            std::cout << "     MOVABLE: ";
        } else {
            std::cout << " NOT MOVABLE: ";
        }
        std::cout << i << " - " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2]
                  << " - " << a.ad << " - " << a.xs << " - " << a.charge << "\n";
    }
}

/**
 * @brief 显示作用在每个原子上的力
 * 
 * 输出每个原子当前受到的力向量（存储在minus_forces中）
 * 
 * @note 用于调试力计算和优化过程
 */
void model::show_forces() const {
    std::cout << "ATOM FORCES\n";
    VINA_FOR_IN(i, atoms) {
        std::cout << i << " " << minus_forces[i][0] << " " << minus_forces[i][1] << " " << minus_forces[i][2] << "\n";
    }
}

/**
 * @brief 综合调试信息输出函数
 * 
 * 根据参数选择性地输出各种调试信息
 * 
 * @param show_coords 是否显示原子坐标
 * @param show_internal 是否显示内部信息（当前未使用）
 * @param show_atoms 是否显示原子详细信息
 * @param show_grid 是否显示网格原子信息
 * @param show_about 是否显示模型统计信息
 * 
 * @note 用于全面的模型调试和分析
 */
void model::print_stuff(bool show_coords, bool show_internal, bool show_atoms, bool show_grid, bool show_about) const {
    
    if (show_coords) {
        std::cout << "coords:\n";
        VINA_FOR_IN(i, coords)
            printnl(coords[i]);  // 输出每个原子的坐标
    }

    if (show_atoms) {
        std::cout << "atoms:\n";
        VINA_FOR_IN(i, atoms) {
            const atom& a = atoms[i];
            // 输出原子的类型信息和电荷
            std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
            std::cout << a.bonds.size() << "  "; printnl(a.coords);  // 输出键数量和坐标
        }
    }

    if (show_grid) {
        std::cout << "grid_atoms:\n";
        VINA_FOR_IN(i, grid_atoms) {
            const atom& a = grid_atoms[i];
            // 输出网格原子的详细信息
            std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
            std::cout << a.bonds.size() << "  "; printnl(a.coords);
        }
    }

    if (show_about) {
        about();  // 调用about()函数显示模型统计信息
    }
}

/**
 * @brief 计算两原子间的成对碰撞惩罚
 * 
 * 基于原子间距离和共价半径计算碰撞惩罚值
 * 
 * @param r 原子间的实际距离
 * @param covalent_r 两原子共价半径之和
 * @return 碰撞惩罚值
 * 
 * @note 惩罚函数特点：
 *       - r = 0时达到最大惩罚
 *       - r = covalent_r时惩罚值为1
 *       - r > 2*covalent_r时惩罚为0
 *       - 其他情况使用双曲函数平滑过渡
 */
fl pairwise_clash_penalty(fl r, fl covalent_r) {
    // r = 0          -> max_penalty 
    // r = covalent_r -> 1
    // elsewhere      -> hyperbolic function
    assert(r >= 0);
    assert(covalent_r > epsilon_fl);
    const fl x = r / covalent_r;
    if(x > 2) return 0;      // 距离足够远时无惩罚
    return 1-x*x/4;          // 双曲衰减函数
}

/**
 * @brief 计算指定相互作用对的碰撞惩罚辅助函数
 * 
 * 遍历给定的相互作用对列表，计算总的碰撞惩罚
 * 
 * @param pairs 要计算的相互作用对列表
 * @return 总的碰撞惩罚值
 * 
 * @note 使用原子的共价半径来定义合理的原子间距离
 */
fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
    fl e = 0;
    VINA_FOR_IN(i, pairs) {
        const interacting_pair& ip = pairs[i];
        const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));  // 计算实际距离
        const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();  // 共价半径和
        e += pairwise_clash_penalty(r, covalent_r);  // 累加碰撞惩罚
    }
    return e;
}

/**
 * @brief 计算模型的总碰撞惩罚
 * 
 * 计算所有配体内部相互作用对和其他相互作用对的碰撞惩罚总和
 * 
 * @return 模型的总碰撞惩罚值
 * 
 * @note 用于评估构象的物理合理性，过高的碰撞惩罚表明构象不合理
 */
fl model::clash_penalty() const {
    fl e = 0;
    VINA_FOR_IN(i, ligands) 
        e += clash_penalty_aux(ligands[i].pairs);  // 累加所有配体的内部碰撞惩罚
    e += clash_penalty_aux(other_pairs);           // 累加其他相互作用对的碰撞惩罚
    return e;
}