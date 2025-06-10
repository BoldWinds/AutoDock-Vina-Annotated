/**
 * @file potentials.h
 * @brief 分子对接势能函数定义
 * 
 * 定义了AutoDock Vina、Vinardo和AD4评分函数中使用的各种势能项，
 * 包括高斯函数、排斥力、疏水相互作用、氢键、静电和溶剂化效应等。

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

#ifndef VINA_POTENTIALS_H
#define VINA_POTENTIALS_H

#include "atom.h"
#include "int_pow.h"

/**
 * @brief 平滑阶跃函数
 * @param x_bad 不利距离
 * @param x_good 有利距离
 * @param x 当前距离
 * @return 0到1之间的平滑过渡值
 * @note 在两个距离值之间提供线性插值，用于势能函数的平滑过渡
 */
inline fl slope_step(fl x_bad, fl x_good, fl x) {
    if (x_bad < x_good) {
        if (x <= x_bad) return 0;
        if (x >= x_good) return 1;
    }
    else {
        if (x >= x_bad) return 0;
        if (x <= x_good) return 1;
    }
    return (x - x_bad) / (x_good - x_bad);
}

/**
 * @brief 检查是否为胶水类型原子
 * @param xs_t XS原子类型
 * @return 如果是G0、G1、G2、G3类型返回true
 * @note 胶水类型用于大环化合物的内部约束
 */
inline bool is_glue_type(sz xs_t) {
    if ((xs_t==XS_TYPE_G0) || (xs_t==XS_TYPE_G1) || (xs_t==XS_TYPE_G2) || (xs_t==XS_TYPE_G3)) return true;
    return false;
}

/**
 * @brief 计算两个原子间的最优距离（Vina标准）
 * @param xs_t1 第一个原子的XS类型
 * @param xs_t2 第二个原子的XS类型
 * @return 最优相互作用距离
 * @note 对于胶水类型返回0，否则返回两个原子半径之和
 */
inline fl optimal_distance(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0; // G0, G1, G2 or G3
    return xs_radius(xs_t1) + xs_radius(xs_t2);
}

/**
 * @brief 安全除法函数
 * @param x 被除数
 * @param y 除数
 * @return 除法结果，处理分母为零的情况
 * @note 避免除零错误，当分母接近零时返回适当的极值
 */
inline fl smooth_div(fl x, fl y) {
    if (std::abs(x) < epsilon_fl) return 0;
    if (std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
    return x / y;
}

/**
 * @brief 计算两个原子间的最优距离（Vinardo标准）
 * @param xs_t1 第一个原子的XS类型
 * @param xs_t2 第二个原子的XS类型
 * @return 最优相互作用距离
 * @note 使用Vinardo特定的原子半径计算
 */
inline fl optimal_distance_vinardo(sz xs_t1, sz xs_t2) {
    if (is_glue_type(xs_t1) || is_glue_type(xs_t2)) return 0.0; // G0, G1, G2 or G3
    return xs_vinardo_radius(xs_t1) + xs_vinardo_radius(xs_t2);
}

/**
 * @brief 距离平滑函数
 * @param r 实际距离
 * @param rij 理想距离
 * @param smoothing 平滑参数
 * @return 平滑后的距离
 * @note 在理想距离附近提供平滑过渡，避免势能函数的不连续性
 */
inline fl smoothen(fl r, fl rij, fl smoothing) {
    fl out;
    smoothing *= 0.5;

    if (r > rij + smoothing)
        out = r - smoothing;
    else if(r < rij - smoothing)
        out = r + smoothing;
    else
        out = rij;

    return out;
}

/**
 * @brief 获取AD4氢键深度参数
 * @param a 原子类型索引
 * @return 氢键深度值
 */
inline fl ad4_hb_eps(sz& a) {
    if (a < AD_TYPE_SIZE) return ad_type_property(a).hb_depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

/**
 * @brief 获取AD4氢键半径参数
 * @param t 原子类型索引
 * @return 氢键半径值
 */
inline fl ad4_hb_radius(sz& t) {
    if (t < AD_TYPE_SIZE) return ad_type_property(t).hb_radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

/**
 * @brief 获取AD4范德华深度参数
 * @param a 原子类型索引
 * @return 范德华深度值
 */
inline fl ad4_vdw_eps(sz& a) {
    if(a < AD_TYPE_SIZE) return ad_type_property(a).depth;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

/**
 * @brief 获取AD4范德华半径参数
 * @param t 原子类型索引
 * @return 范德华半径值
 */
inline fl ad4_vdw_radius(sz& t) {
    if(t < AD_TYPE_SIZE) return ad_type_property(t).radius;
    VINA_CHECK(false);
    return 0; // placating the compiler
}

/**
 * @brief 检查两个原子是否通过胶水类型连接
 * @param xs_t1 第一个原子的XS类型
 * @param xs_t2 第二个原子的XS类型
 * @return 如果两原子通过胶水类型连接返回true
 * @note 用于大环化合物中特定原子类型对的识别
 */
inline bool is_glued(sz xs_t1, sz xs_t2) {
    return (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_H_CG0) ||
       (xs_t1 == XS_TYPE_G0 && xs_t2 == XS_TYPE_C_P_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_H_CG0) ||
       (xs_t2 == XS_TYPE_G0 && xs_t1 == XS_TYPE_C_P_CG0) ||

       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_H_CG1) ||
       (xs_t1 == XS_TYPE_G1 && xs_t2 == XS_TYPE_C_P_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_H_CG1) ||
       (xs_t2 == XS_TYPE_G1 && xs_t1 == XS_TYPE_C_P_CG1) ||

       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_H_CG2) ||
       (xs_t1 == XS_TYPE_G2 && xs_t2 == XS_TYPE_C_P_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_H_CG2) ||
       (xs_t2 == XS_TYPE_G2 && xs_t1 == XS_TYPE_C_P_CG2) ||

       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_H_CG3) ||
       (xs_t1 == XS_TYPE_G3 && xs_t2 == XS_TYPE_C_P_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_H_CG3) ||
       (xs_t2 == XS_TYPE_G3 && xs_t1 == XS_TYPE_C_P_CG3);
}

/**
 * @class Potential
 * @brief 势能函数基类
 * @note 为所有具体势能函数提供统一接口
 */
class Potential {
public:
    virtual ~Potential() { }
    /**
     * @brief 计算两个原子间的势能
     * @param a 第一个原子
     * @param b 第二个原子
     * @param r 原子间距离
     * @return 势能值
     */
    virtual fl eval(const atom& a, const atom& b, fl r) { return 0; };
    
    /**
     * @brief 根据原子类型计算势能
     * @param t1 第一个原子类型
     * @param t2 第二个原子类型
     * @param r 原子间距离
     * @return 势能值
     */
    virtual fl eval(sz t1, sz t2, fl r) { return 0; };
    
    /**
     * @brief 获取势能函数的截断距离
     * @return 截断距离
     */
    virtual fl get_cutoff() { return 0; }
};

/**
 * @class vina_gaussian
 * @brief Vina高斯势能函数
 * @note 用于模拟有利的分子间相互作用，在最优距离处达到最大值
 */
class vina_gaussian : public Potential {
public:
    /**
     * @brief 构造函数
     * @param offset_ 距离偏移量
     * @param width_ 高斯函数宽度
     * @param cutoff_ 截断距离
     */
    vina_gaussian(fl offset_, fl width_, fl cutoff_) : offset(offset_), width(width_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        return gauss(r - (optimal_distance(a.xs, b.xs) + offset)); // hard-coded to XS atom type
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        return gauss(r - (optimal_distance(t1, t2) + offset)); // hard-coded to XS atom type
    };

    fl get_cutoff() { return cutoff; }

private:
    fl offset; ///< 添加到最优距离的偏移量
    fl width;  ///< 高斯函数的宽度参数
    fl cutoff; ///< 截断距离

    /**
     * @brief 高斯函数计算
     * @param x 输入值
     * @return 高斯函数值
     * @note 计算 exp(-(x/width)^2)
     */
    fl gauss(fl x) {
        return std::exp(-sqr(x / width));
    };
};

/**
 * @class vina_repulsion
 * @brief Vina排斥势能函数
 * @note 当原子距离小于最优距离时产生排斥力，防止原子重叠
 */
class vina_repulsion : public Potential {
public:
    /**
     * @brief 构造函数
     * @param offset_ 距离偏移量
     * @param cutoff_ 截断距离
     */
    vina_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        fl d = r - (optimal_distance(a.xs, b.xs) + offset); // hard-coded to XS atom type
        if (d > 0.0)
            return 0.0;
        return d * d; // 距离小于最优距离时产生二次排斥
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl d = r - (optimal_distance(t1, t2) + offset); // hard-coded to XS atom type
        if(d > 0.0)
            return 0.0;
        return d*d;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl offset; ///< 添加到范德华半径的偏移量
    fl cutoff; ///< 截断距离
};

/**
 * @class vina_hydrophobic
 * @brief Vina疏水相互作用势能函数
 * @note 当两个疏水原子在适当距离时产生有利相互作用
 */
class vina_hydrophobic : public Potential {
public:
    /**
     * @brief 构造函数
     * @param good_ 有利距离
     * @param bad_ 不利距离
     * @param cutoff_ 截断距离
     */
    vina_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        else return 0.0;
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        else
            return 0.0;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl good;   ///< 有利相互作用距离
    fl bad;    ///< 不利相互作用距离
    fl cutoff; ///< 截断距离
};

/**
 * @class vina_non_dir_h_bond
 * @brief Vina非定向氢键势能函数
 * @note 计算两个可能形成氢键的原子间的相互作用
 */
class vina_non_dir_h_bond : public Potential {
public:
    /**
     * @brief 构造函数
     * @param good_ 有利距离
     * @param bad_ 不利距离
     * @param cutoff_ 截断距离
     */
    vina_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs))
            return slope_step(bad, good, r - optimal_distance(a.xs, b.xs));
        return 0.0;
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance(t1, t2));
        return 0.0;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl good;   ///< 有利氢键距离
    fl bad;    ///< 不利氢键距离
    fl cutoff; ///< 截断距离
};

// Vinardo
/**
 * @class vinardo_gaussian
 * @brief Vinardo高斯势能函数
 * @note 用于模拟有利的分子间相互作用，在最优距离处达到最大值
 */
class vinardo_gaussian : public Potential {
public:
    /**
     * @brief 构造函数
     * @param offset_ 距离偏移量
     * @param width_ 高斯函数宽度
     * @param cutoff_ 截断距离
     */
    vinardo_gaussian(fl offset_, fl width_, fl cutoff_) : offset(offset_), width(width_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        return gauss(r - (optimal_distance_vinardo(a.xs, b.xs) + offset)); // hard-coded to XS atom type
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        return gauss(r - (optimal_distance_vinardo(t1, t2) + offset)); // hard-coded to XS atom type
    };

    fl get_cutoff() { return cutoff; }

private:
    fl offset; ///< 添加到最优距离的偏移量
    fl width;  ///< 高斯函数的宽度参数
    fl cutoff; ///< 截断距离

    /**
     * @brief 高斯函数计算
     * @param x 输入值
     * @return 高斯函数值
     * @note 计算 exp(-(x/width)^2)
     */
    fl gauss(fl x) {
        return std::exp(-sqr(x / width));
    };
};

/**
 * @class vinardo_repulsion
 * @brief Vinardo排斥势能函数
 * @note 当原子距离小于最优距离时产生排斥力，防止原子重叠
 */
class vinardo_repulsion : public Potential {
public:
    /**
     * @brief 构造函数
     * @param offset_ 距离偏移量
     * @param cutoff_ 截断距离
     */
    vinardo_repulsion(fl offset_, fl cutoff_) : offset(offset_), cutoff(cutoff_) {}

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        fl d = r - (optimal_distance_vinardo(a.xs, b.xs) + offset); // hard-coded to XS atom type
        if (d > 0.0)
            return 0.0;
        return d * d; // 距离小于最优距离时产生二次排斥
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl d = r - (optimal_distance_vinardo(t1, t2) + offset); // hard-coded to XS atom type
        if(d > 0.0)
            return 0.0;
        return d*d;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl offset; ///< 添加到范德华半径的偏移量
    fl cutoff; ///< 截断距离
};

/**
 * @class vinardo_hydrophobic
 * @brief Vinardo疏水相互作用势能函数
 * @note 当两个疏水原子在适当距离时产生有利相互作用
 */
class vinardo_hydrophobic : public Potential {
public:
    /**
     * @brief 构造函数
     * @param good_ 有利距离
     * @param bad_ 不利距离
     * @param cutoff_ 截断距离
     */
    vinardo_hydrophobic(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_is_hydrophobic(a.xs) && xs_is_hydrophobic(b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        else return 0.0;
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        else
            return 0.0;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl good;   ///< 有利相互作用距离
    fl bad;    ///< 不利相互作用距离
    fl cutoff; ///< 截断距离
};

/**
 * @class vinardo_non_dir_h_bond
 * @brief Vinardo非定向氢键势能函数
 * @note 计算两个可能形成氢键的原子间的相互作用
 */
class vinardo_non_dir_h_bond : public Potential {
public:
    /**
     * @brief 构造函数
     * @param good_ 有利距离
     * @param bad_ 不利距离
     * @param cutoff_ 截断距离
     */
    vinardo_non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : good(good_), bad(bad_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if ((a.xs >= XS_TYPE_SIZE) || (b.xs >= XS_TYPE_SIZE))
            return 0.0;
        if (xs_h_bond_possible(a.xs, b.xs))
            return slope_step(bad, good, r - optimal_distance_vinardo(a.xs, b.xs));
        return 0.0;
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if(xs_h_bond_possible(t1, t2))
            return slope_step(bad, good, r - optimal_distance_vinardo(t1, t2));
        return 0.0;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl good;   ///< 有利氢键距离
    fl bad;    ///< 不利氢键距离
    fl cutoff; ///< 截断距离
};

// AD42
/**
 * @class ad4_electrostatic
 * @brief AD4静电势能函数
 * @note 计算带电原子间的静电相互作用
 */
class ad4_electrostatic : public Potential {
public:
    ad4_electrostatic(fl cap_, fl cutoff_) : cap(cap_), cutoff(cutoff_) { }
    //~ad4_electrostatic() { }
    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl q1q2 = a.charge * b.charge * 332.0;
        fl B = 78.4 + 8.5525;
        fl lB = -B * 0.003627;
        fl diel = -8.5525 + (B / (1 + 7.7839 * std::exp(lB * r)));
        if (r < epsilon_fl)
            return q1q2 * cap / diel;
        else {
            return q1q2 * (std::min)(cap, 1.0 / (r * diel));
        }
    };
    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }
private:
    fl cap;
    fl cutoff;
};

/**
 * @class ad4_solvation
 * @brief AD4溶剂化势能函数
 * @note 计算原子的溶剂化自由能贡献，基于原子体积和溶剂化参数
 */
class ad4_solvation : public Potential {
public:
    /**
     * @brief 构造函数
     * @param desolvation_sigma_ 去溶剂化高斯宽度
     * @param solvation_q_ 电荷相关溶剂化参数
     * @param charge_dependent_ 是否考虑电荷依赖性
     * @param cutoff_ 截断距离
     */
    ad4_solvation(fl desolvation_sigma_, fl solvation_q_, bool charge_dependent_, fl cutoff_) : solvation_q(solvation_q_), charge_dependent(charge_dependent_), desolvation_sigma(desolvation_sigma_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        fl q1 = a.charge;
        fl q2 = b.charge;
        VINA_CHECK(not_max(q1));
        VINA_CHECK(not_max(q2));
        fl solv1 = solvation_parameter(a);
        fl solv2 = solvation_parameter(b);
        fl volume1 = volume(a);
        fl volume2 = volume(b);
        fl my_solv = charge_dependent ? solvation_q : 0;
        // 计算溶剂化能：考虑原子体积、溶剂化参数和距离衰减
        fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 + 
                  (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-0.5 * sqr(r / desolvation_sigma));
        VINA_CHECK(not_max(tmp));
        return tmp;
    };

    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }

private:
    fl desolvation_sigma; ///< 去溶剂化高斯函数的标准差
    fl solvation_q;       ///< 电荷相关的溶剂化参数
    bool charge_dependent; ///< 是否考虑电荷依赖的溶剂化
    fl cutoff;            ///< 截断距离

    /**
     * @brief 计算原子体积
     * @param a 原子类型
     * @return 原子体积
     * @note 优先使用AD类型的体积参数，否则基于XS半径计算球体体积
     */
    fl volume(const atom_type& a) const {
        if (a.ad < AD_TYPE_SIZE)
            return ad_type_property(a.ad).volume;
        else if (a.xs < XS_TYPE_SIZE)
            return 4.0 * pi / 3.0 * int_pow<3>(xs_radius(a.xs));
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };

    /**
     * @brief 获取原子的溶剂化参数
     * @param a 原子类型
     * @return 溶剂化参数
     * @note 对金属原子使用特殊的溶剂化参数
     */
    fl solvation_parameter(const atom_type& a) const {
        if (a.ad < AD_TYPE_SIZE)
            return ad_type_property(a.ad).solvation;
        else if (a.xs == XS_TYPE_Met_D)
            return metal_solvation_parameter;
        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };
};

/**
 * @class ad4_vdw
 * @brief AD4范德华势能函数
 * @note 实现12-6 Lennard-Jones势，用于计算范德华相互作用
 */
class ad4_vdw : public Potential {
public:
    /**
     * @brief 构造函数
     * @param smoothing_ 平滑参数
     * @param cap_ 势能上限
     * @param cutoff_ 截断距离
     */
    ad4_vdw(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);
        fl vdw_depth = std::sqrt(ad4_vdw_eps(t1) * ad4_vdw_eps(t2));

        if (hb_depth < 0) return 0.0; // 如果是氢键相互作用，则不计算范德华力

        r = smoothen(r, vdw_rij, smoothing); // 平滑距离

        // 计算12-6 Lennard-Jones势的系数
        fl c_12 = int_pow<12>(vdw_rij) * vdw_depth;
        fl c_6  = int_pow<6>(vdw_rij)  * vdw_depth * 2.0;
        fl r6   = int_pow<6>(r);
        fl r12  = int_pow<12>(r);

        if(r12 > epsilon_fl && r6 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_6 / r6); // U = C12/r^12 - C6/r^6
        else
            return cap; // 避免除零，返回上限值

        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };

    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }

private:
    fl smoothing; ///< 距离平滑参数
    fl cap;       ///< 势能上限
    fl cutoff;    ///< 截断距离
};

/**
 * @class ad4_hb
 * @brief AD4氢键势能函数
 * @note 实现12-10势函数来模拟氢键相互作用
 */
class ad4_hb : public Potential {
public:
    /**
     * @brief 构造函数
     * @param smoothing_ 平滑参数
     * @param cap_ 势能上限
     * @param cutoff_ 截断距离
     */
    ad4_hb(fl smoothing_, fl cap_, fl cutoff_) : smoothing(smoothing_), cap(cap_), cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        sz t1 = a.ad;
        sz t2 = b.ad;
        fl hb_rij = ad4_hb_radius(t1) + ad4_hb_radius(t2);
        fl hb_depth = ad4_hb_eps(t1) * ad4_hb_eps(t2);
        fl vdw_rij = ad4_vdw_radius(t1) + ad4_vdw_radius(t2);

        if (hb_depth >= 0)
            return 0.0; // 如果不是氢键相互作用（深度非负），则返回0

        r = smoothen(r, hb_rij, smoothing); // 平滑距离

        // 计算12-10氢键势的系数
        fl c_12 = int_pow<12>(hb_rij) * -hb_depth * 10 / 2.0;
        fl c_10 = int_pow<10>(hb_rij) * -hb_depth * 12 / 2.0;
        fl r10  = int_pow<10>(r);
        fl r12  = int_pow<12>(r);

        if (r12 > epsilon_fl && r10 > epsilon_fl)
            return (std::min)(cap, c_12 / r12 - c_10 / r10); // U = C12/r^12 - C10/r^10
        else
            return cap; // 避免除零，返回上限值

        VINA_CHECK(false);
        return 0.0; // placating the compiler
    };

    fl eval(sz t1, sz t2, fl r) { return 0; }
    fl get_cutoff() { return cutoff; }

private:
    fl smoothing; ///< 距离平滑参数
    fl cap;       ///< 势能上限
    fl cutoff;    ///< 截断距离
};

/**
 * @class linearattraction
 * @brief 线性吸引势能函数
 * @note 用于大环化合物中胶水原子对的线性约束，势能与距离成正比
 */
class linearattraction : public Potential {
public:
    /**
     * @brief 构造函数
     * @param cutoff_ 截断距离
     */
    linearattraction(fl cutoff_): cutoff(cutoff_) { }

    fl eval(const atom& a, const atom& b, fl r) {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(a.xs, b.xs))
            return r; // 对于胶水原子对，势能与距离成正比
        else 
            return 0.0;
    };

    fl eval(sz t1, sz t2, fl r) {
        if (r >= cutoff)
            return 0.0;
        if (is_glued(t1, t2))
            return r;
        else
            return 0.0;
    };

    fl get_cutoff() { return cutoff; }

private:
    fl cutoff; ///< 截断距离
};

#endif
