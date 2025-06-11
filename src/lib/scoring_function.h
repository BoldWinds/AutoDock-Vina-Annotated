/**
 * @file scoring_function.h
 * @brief AutoDock Vina评分函数系统实现
 * 
 * 该文件实现了分子对接中的能量评分计算功能，支持多种评分函数类型。

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

#ifndef VINA_SCORING_FUNCTION_H
#define VINA_SCORING_FUNCTION_H

#include <stdlib.h>
#include <list>
#include "atom.h"
#include "conf_independent.h"
#include "potentials.h"
#include "common.h"

//Forward declaration
struct model;

/**
 * @brief 评分函数类型枚举
 * 
 * 定义了AutoDock Vina支持的三种不同评分函数
 */
enum scoring_function_choice {
    SF_VINA,    ///< Vina原始评分函数
    SF_AD42,    ///< AutoDock 4.2评分函数
    SF_VINARDO  ///< Vinardo评分函数
};

/**
 * @brief 分子对接评分函数类
 * 
 * 该类实现了分子对接中的能量评分计算，包含多个势函数和构象无关项。
 * 根据不同的评分函数类型配置相应的势函数组合和参数。
 */
class ScoringFunction {
public:
    /**
     * @brief 默认构造函数
     * 
     * 初始化势函数和构象无关项数量为0
     */
    ScoringFunction() {
        m_num_potentials = 0;
        m_num_conf_independents = 0;
    }

    /**
     * @brief 带参数的构造函数
     * 
     * 根据指定的评分函数类型和权重初始化评分函数
     * 
     * @param sf_choice 评分函数类型选择
     * @param weights 各势函数的权重向量
     * 
     * @note 不同评分函数使用不同的势函数组合：
     *       - VINA: 高斯势、排斥势、疏水作用、氢键、线性吸引势
     *       - VINARDO: 改进的高斯势和其他势函数参数
     *       - AD42: AutoDock 4.2的vdW、氢键、静电和溶剂化势
     */
    ScoringFunction(const scoring_function_choice sf_choice, const flv& weights){
        switch (sf_choice)
        {
            case SF_VINA:
            {
                // Vina评分函数：添加各种势函数
                m_potentials.push_back(new vina_gaussian(0, 0.5, 8.0));        ///< 高斯势1（offset 0Å，width 0.5Å）
                m_potentials.push_back(new vina_gaussian(3, 2.0, 8.0));        ///< 高斯势2（offset 3Å，width 2.0Å）
                m_potentials.push_back(new vina_repulsion(0.0, 8.0));          ///< 排斥势
                m_potentials.push_back(new vina_hydrophobic(0.5, 1.5, 8.0));   ///< 疏水作用势
                m_potentials.push_back(new vina_non_dir_h_bond(-0.7, 0, 8.0)); ///< 非定向氢键势
                m_potentials.push_back(new linearattraction(20.0));            ///< 线性吸引势
                m_conf_independents.push_back(new num_tors_div());             ///< 扭转角惩罚项
                m_atom_typing = atom_type::XS;  ///< 使用XS原子类型
                m_cutoff = 8.0;                 ///< 截断距离8Å
                m_max_cutoff = 20.0;            ///< 最大截断距离20Å
                break;
            }
            case SF_VINARDO:
            {
                // VINARDO评分函数：使用改进的势函数参数
                m_potentials.push_back(new vinardo_gaussian(0, 0.8, 8.0));        ///< VINARDO高斯势
                m_potentials.push_back(new vinardo_repulsion(0, 8.0));            ///< VINARDO排斥势
                m_potentials.push_back(new vinardo_hydrophobic(0, 2.5, 8.0));     ///< VINARDO疏水势
                m_potentials.push_back(new vinardo_non_dir_h_bond(-0.6, 0, 8.0)); ///< VINARDO氢键势
                m_potentials.push_back(new linearattraction(20.0));               ///< 线性吸引势
                m_conf_independents.push_back(new num_tors_div());                ///< 扭转角惩罚项
                m_atom_typing = atom_type::XS;
                m_cutoff = 8.0;
                m_max_cutoff = 20.0;
                break;
            }
            case SF_AD42:
            {
                // AutoDock 4.2评分函数：使用AD4势函数
                m_potentials.push_back(new ad4_vdw(0.5, 100000, 8.0));           ///< AD4范德华势
                m_potentials.push_back(new ad4_hb(0.5, 100000, 8.0));            ///< AD4氢键势
                m_potentials.push_back(new ad4_electrostatic(100, 20.48));        ///< AD4静电势
                m_potentials.push_back(new ad4_solvation(3.6, 0.01097, true, 20.48)); ///< AD4溶剂化势
                m_potentials.push_back(new linearattraction(20.0));               ///< 线性吸引势
                m_conf_independents.push_back(new ad4_tors_add());                ///< AD4扭转角项
                m_atom_typing = atom_type::AD;  ///< 使用AD原子类型
                m_cutoff = 20.48;               ///< 截断距离20.48Å
                m_max_cutoff = 20.48;
                break;
            }
            default:
            {
                std::cout << "INSIDE everything::everything()   sfchoice = " << sf_choice << "\n";
                VINA_CHECK(false); ///< 未知评分函数类型，触发错误
                break;
            }
        }

        // 设置势函数和构象无关项的数量
        m_num_potentials = m_potentials.size();
        m_num_conf_independents = m_conf_independents.size();
        m_weights = weights; ///< 保存权重向量
    };

    /**
     * @brief 销毁所有动态分配的势函数和构象无关项
     * 
     * 释放内存并重置计数器
     */
    void Destroy()
    {
        // 删除所有势函数对象
        for (auto p : m_potentials)
        {
            delete p;
        }
        m_potentials.clear();
        m_num_potentials = 0;

        // 删除所有构象无关项对象
        for (auto p : m_conf_independents)
        {
            delete p;
        }
        m_conf_independents.clear();
        m_num_conf_independents = 0;
    }

    /**
     * @brief 析构函数
     * 
     * 调用Destroy()释放所有资源
     */
    ~ScoringFunction() {
        Destroy();
    }

    /**
     * @brief 计算两个原子间的相互作用能
     * 
     * @param a 第一个原子
     * @param b 第二个原子
     * @param r 原子间距离
     * @return fl 加权后的总相互作用能
     * 
     * @note 该函数不检查截断距离，调用者需要确保距离在有效范围内
     */
    fl eval(atom& a, atom& b, fl r) const{
        fl acc = 0;
        // 遍历所有势函数，计算加权能量总和
        VINA_FOR (i, m_num_potentials)
        {
            acc += m_weights[i] * m_potentials[i]->eval(a, b, r);
        }
        return acc;
    };

    /**
     * @brief 根据原子类型索引计算相互作用能
     * 
     * @param t1 第一个原子的类型索引
     * @param t2 第二个原子的类型索引
     * @param r 原子间距离
     * @return fl 加权后的总相互作用能
     * @note 似乎是在预计算的地方发挥作用？
     */
    fl eval(sz t1, sz t2, fl r) const{
        fl acc = 0;
        // 遍历所有势函数，计算加权能量总和
        VINA_FOR (i, m_num_potentials)
        {
            acc += m_weights[i] * m_potentials[i]->eval(t1, t2, r);
        }
        return acc;
    };

    /**
     * @brief 计算构象无关项的贡献
     * 
     * @param m 分子模型
     * @param e 输入的能量值
     * @return fl 经过构象无关项修正后的能量
     * 
     * @note 构象无关项通常包括扭转角惩罚等，不直接累加而是对能量进行修正
     */
    fl conf_independent(const model& m, fl e) const{
        // 权重迭代器，指向势函数权重之后的构象无关项权重
        flv::const_iterator it = m_weights.begin() + m_num_potentials;
        conf_independent_inputs in(m); // 构造输入参数（效率较低但此处速度不重要）
        
        VINA_FOR (i, m_num_conf_independents)
        {
            // 注意：这里不累加能量，而是对能量进行修正
            e = m_conf_independents[i]->eval(in, e, it);
        }
        assert(it == m_weights.end()); // 确保所有权重都被使用
        return e;
    };

    /**
     * @brief 获取截断距离
     * @return fl 截断距离值
     */
    fl get_cutoff() const { return m_cutoff; }

    /**
     * @brief 获取最大截断距离
     * @return fl 最大截断距离值
     */
    fl get_max_cutoff() const { return m_max_cutoff; }

    /**
     * @brief 获取原子类型系统
     * @return atom_type::t 原子类型枚举值
     */
    atom_type::t get_atom_typing() const { return m_atom_typing; }

    /**
     * @brief 获取所有原子类型的索引向量
     * @return szv 包含所有原子类型索引的向量
     */
    szv get_atom_types() const{
        szv tmp;
        // 遍历所有原子类型，添加到向量中
        VINA_FOR(i, num_atom_types(m_atom_typing))
        {
          tmp.push_back(i);
        }
        return tmp;
    };

    /**
     * @brief 获取原子类型的总数
     * @return sz 原子类型数量
     */
    sz get_num_atom_types() const { return num_atom_types(m_atom_typing); }

    /**
     * @brief 获取权重向量
     * @return flv 包含所有势函数和构象无关项权重的向量
     */
    flv get_weights() const { return m_weights; }

private:
    std::vector<Potential*> m_potentials;           ///< 势函数指针向量
    std::vector<ConfIndependent*> m_conf_independents; ///< 构象无关项指针向量
    flv m_weights;                                  ///< 权重向量
    fl m_cutoff;                                    ///< 截断距离
    fl m_max_cutoff;                                ///< 最大截断距离
    int m_num_potentials;                           ///< 势函数数量
    int m_num_conf_independents;                    ///< 构象无关项数量
    atom_type::t m_atom_typing;                     ///< 原子类型系统
};

#endif
