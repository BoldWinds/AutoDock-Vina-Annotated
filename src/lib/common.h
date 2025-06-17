/**
 * @file common.h
 * @brief Vina常用类型和函数定义

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

#ifndef VINA_COMMON_H
#define VINA_COMMON_H

#include <cassert>
#include <string>
#include <limits>
#include <utility> // pair
#include <algorithm> // too common
#include <vector> // used in typedef, and commonly used overall
#include <cmath> // commonly used
#include <iostream> // various debugging everywhere
#include <fstream> // print_coords
#include <iomanip> // to_string
#include <sstream> // to_string
#include <string> // probably included by the above anyway, common anyway

#include <boost/serialization/vector.hpp> // can't come before the above two - wart fixed in upcoming Boost versions
#include <boost/serialization/base_object.hpp> // movable_atom needs it - (derived from atom)
#include <boost/filesystem/path.hpp> // typedef'ed

#include "macros.h"

typedef double fl;

/**
 * @brief 计算数值的平方
 * @tparam T 数值类型
 * @param x 输入数值
 * @return 返回x的平方
 */
template<typename T>
T sqr(T x) {
    return x*x;
}

/** @brief 非数值常量，用于表示无效或未初始化的数值 */
const fl not_a_num = std::sqrt(fl(-1)); // FIXME? check 

/** @brief 大小类型定义 */
typedef std::size_t sz;

/** @brief 浮点数对类型定义 */
typedef std::pair<fl, fl> pr;

/**
 * @brief 三维向量结构体
 * 
 * 表示三维空间中的向量，提供基本的向量运算功能
 */
struct vec {
    fl data[3];  ///< 存储x,y,z三个分量的数组
    
    /**
     * @brief 默认构造函数
     * 在调试模式下初始化为非数值
     */
    vec() {
#ifndef NDEBUG
        data[0] = data[1] = data[2] = not_a_num;
#endif
    }
    
    /**
     * @brief 参数构造函数
     * @param x X分量
     * @param y Y分量  
     * @param z Z分量
     */
    vec(fl x, fl y, fl z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    
    /**
     * @brief 常量下标访问操作符
     * @param i 分量索引(0,1,2)
     * @return 返回第i个分量的常量引用
     */
    const fl& operator[](sz i) const { assert(i < 3); return data[i]; }
    
    /**
     * @brief 非常量下标访问操作符
     * @param i 分量索引(0,1,2)
     * @return 返回第i个分量的引用
     */
    fl& operator[](sz i) { assert(i < 3); return data[i]; }
    
    /**
     * @brief 计算向量长度的平方
     * @return 返回|v|²
     */
    fl norm_sqr() const {
        return sqr(data[0]) + sqr(data[1]) + sqr(data[2]);
    }
    
    /**
     * @brief 计算向量的长度(模)
     * @return 返回|v|
     */
    fl norm() const {
        return std::sqrt(norm_sqr());
    }
    
    /**
     * @brief 向量点积运算
     * @param v 另一个向量
     * @return 返回两向量的点积
     */
    fl operator*(const vec& v) const {
        return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
    }
    
    /**
     * @brief 向量加法赋值运算
     * @param v 要加的向量
     * @return 返回自身引用
     */
    const vec& operator+=(const vec& v) {
        data[0] += v[0];
        data[1] += v[1];
        data[2] += v[2];
        return *this;
    }
    
    /**
     * @brief 向量减法赋值运算
     * @param v 要减的向量
     * @return 返回自身引用
     */
    const vec& operator-=(const vec& v) {
        data[0] -= v[0];
        data[1] -= v[1];
        data[2] -= v[2];
        return *this;
    }
    
    /**
     * @brief 标量加法赋值运算(每个分量都加上标量)
     * @param s 标量值
     * @return 返回自身引用
     */
    const vec& operator+=(fl s) {
        data[0] += s;
        data[1] += s;
        data[2] += s;
        return *this;
    }
    
    /**
     * @brief 标量减法赋值运算(每个分量都减去标量)
     * @param s 标量值
     * @return 返回自身引用
     */
    const vec& operator-=(fl s) {
        data[0] -= s;
        data[1] -= s;
        data[2] -= s;
        return *this;
    }
    
    /**
     * @brief 向量加法运算
     * @param v 要加的向量
     * @return 返回新的向量结果
     */
    vec operator+(const vec& v) const {
        return vec(data[0] + v[0],
                   data[1] + v[1],
                   data[2] + v[2]);
    }
    
    /**
     * @brief 向量减法运算
     * @param v 要减的向量
     * @return 返回新的向量结果
     */
    vec operator-(const vec& v) const {
        return vec(data[0] - v[0],
                   data[1] - v[1],
                   data[2] - v[2]);
    }
    
    /**
     * @brief 标量乘法赋值运算
     * @param s 标量值
     * @return 返回自身引用
     */
    const vec& operator*=(fl s) {
        data[0] *= s;
        data[1] *= s;
        data[2] *= s;
        return *this;
    }
    
    /**
     * @brief 将所有分量设置为指定值
     * @param s 要设置的值
     */
    void assign(fl s) {
        data[0] = data[1] = data[2] = s;
    }
    
    /**
     * @brief 返回向量维度
     * @return 始终返回3
     */
    sz size() const { return 3; }
    
private:
    friend class boost::serialization::access;
    template<class Archive> 
    void serialize(Archive& ar, const unsigned version) {
        ar & data;
    }
};

/**
 * @brief 标量与向量的乘法运算
 * @param s 标量值
 * @param v 向量
 * @return 返回标量乘以向量的结果
 */
inline vec operator*(fl s, const vec& v) {
    return vec(s * v[0], s * v[1], s * v[2]);
}

/**
 * @brief 计算两个向量的叉积
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 返回a×b的叉积结果
 */
inline vec cross_product(const vec& a, const vec& b) {
    return vec(a[1]*b[2] - a[2]*b[1],
               a[2]*b[0] - a[0]*b[2],
               a[0]*b[1] - a[1]*b[0]);
}

/**
 * @brief 计算两个向量的逐元素乘积
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 返回(a[0]*b[0], a[1]*b[1], a[2]*b[2])
 */
inline vec elementwise_product(const vec& a, const vec& b) {
    return vec(a[0] * b[0],
               a[1] * b[1],
               a[2] * b[2]);
}

/**
 * @brief 3x3矩阵结构体
 * 
 * 使用列主序存储的3x3矩阵，主要用于旋转变换
 */
struct mat {
    fl data[9];  ///< 矩阵数据，列主序存储
    
    /**
     * @brief 默认构造函数
     * 在调试模式下初始化为非数值
     */
    mat() {
#ifndef NDEBUG
        data[0] = data[1] = data[2] =
        data[3] = data[4] = data[5] =
        data[6] = data[7] = data[8] = not_a_num;
#endif
    }
    
    /**
     * @brief 矩阵元素访问(常量版本)
     * @param i 行索引
     * @param j 列索引
     * @return 返回矩阵元素的常量引用
     */
    const fl& operator()(sz i, sz j) const { assert(i < 3); assert(j < 3); return data[i + 3*j]; }
    
    /**
     * @brief 矩阵元素访问(非常量版本)
     * @param i 行索引
     * @param j 列索引
     * @return 返回矩阵元素的引用
     */
    fl& operator()(sz i, sz j) { assert(i < 3); assert(j < 3); return data[i + 3*j]; }

    /**
     * @brief 参数构造函数
     * @param xx,xy,xz 第一行元素
     * @param yx,yy,yz 第二行元素
     * @param zx,zy,zz 第三行元素
     */
    mat(fl xx, fl xy, fl xz,
        fl yx, fl yy, fl yz,
        fl zx, fl zy, fl zz) {

        data[0] = xx; data[3] = xy; data[6] = xz;
        data[1] = yx; data[4] = yy; data[7] = yz;
        data[2] = zx; data[5] = zy; data[8] = zz;
    }
    
    /**
     * @brief 矩阵与向量的乘法运算
     * @param v 输入向量
     * @return 返回矩阵乘以向量的结果
     */
    vec operator*(const vec& v) const {
        return vec(data[0]*v[0] + data[3]*v[1] + data[6]*v[2], 
                   data[1]*v[0] + data[4]*v[1] + data[7]*v[2],
                   data[2]*v[0] + data[5]*v[1] + data[8]*v[2]);
    }
    
    /**
     * @brief 标量乘法赋值运算
     * @param s 标量值
     * @return 返回自身引用
     */
    const mat& operator*=(fl s) {
        VINA_FOR(i, 9)
            data[i] *= s;
        return *this;
    }
};

/** @brief 向量容器类型定义 */
typedef std::vector<vec> vecv;

/** @brief 向量对类型定义 */
typedef std::pair<vec, vec> vecp;

/** @brief 浮点数容器类型定义 */
typedef std::vector<fl> flv;

/** @brief 浮点数对容器类型定义 */
typedef std::vector<pr> prv;

/** @brief 大小类型容器定义 */
typedef std::vector<sz> szv;

/** @brief 文件路径类型定义 */
typedef boost::filesystem::path path;

/**
 * @brief 内部错误异常结构体
 * 
 * 用于记录发生错误的文件和行号
 */
struct internal_error {
    std::string file;  ///< 发生错误的文件名
    unsigned line;     ///< 发生错误的行号
    
    /**
     * @brief 构造函数
     * @param file_ 文件名
     * @param line_ 行号
     */
    internal_error(const std::string& file_, unsigned line_) : file(file_), line(line_) {}
};

#ifdef NDEBUG
    /** @brief 发布模式下的检查宏，失败时抛出异常 */
    #define VINA_CHECK(P) do { if(!(P)) throw internal_error(__FILE__, __LINE__); } while(false)
#else
    /** @brief 调试模式下的检查宏，失败时断言 */
    #define VINA_CHECK(P) assert(P)
#endif

/** @brief 圆周率常量 */
const fl pi = fl(3.1415926535897931);

/**
 * @brief 将浮点数转换为大小类型，限制在指定范围内
 * @param x 输入浮点数
 * @param max_sz 最大值限制
 * @return 返回[0, max_sz]范围内的大小值
 */
inline sz fl_to_sz(fl x, sz max_sz) {
    if(x <= 0) return 0;
    if(x >= max_sz) return max_sz;
    sz tmp = static_cast<sz>(x);
    return (std::min)(tmp, max_sz);
}

/** @brief 浮点数比较容差 */
const fl fl_tolerance = fl(0.001);

/**
 * @brief 判断两个浮点数是否相等(在容差范围内)
 * @param a 第一个数
 * @param b 第二个数
 * @return 如果相等返回true
 */
inline bool eq(fl a, fl b) {
    return std::abs(a - b) < fl_tolerance; 
}

/**
 * @brief 判断两个向量是否相等(在容差范围内)
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 如果相等返回true
 */
inline bool eq(const vec& a, const vec& b) {
    return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

/**
 * @brief 判断两个容器是否相等
 * @tparam T 容器元素类型
 * @param a 第一个容器
 * @param b 第二个容器
 * @return 如果相等返回true
 */
template<typename T>
bool eq(const std::vector<T>& a, const std::vector<T>& b) {
    if(a.size() != b.size()) return false;
    VINA_FOR_IN(i, a)
        if(!eq(a[i], b[i]))
            return false;
    return true;
}

/** @brief 浮点数最大值 */
const fl max_fl = (std::numeric_limits<fl>::max)();

/** @brief 大小类型最大值 */
const sz max_sz = (std::numeric_limits<sz>::max)();

/** @brief 无符号整数最大值 */
const unsigned max_unsigned = (std::numeric_limits<unsigned>::max)();

/** @brief 浮点数机器精度 */
const fl epsilon_fl = std::numeric_limits<fl>::epsilon();

/** @brief 零向量常量 */
const vec zero_vec(0, 0, 0);

/** @brief 最大值向量常量 */
const vec max_vec(max_fl, max_fl, max_fl);

/**
 * @brief 检查浮点数是否不是最大值
 * @param x 待检查的浮点数
 * @return 如果不是最大值返回true
 */
inline bool not_max(fl x) {
    return (x < 0.1 * max_fl);
}

/**
 * @brief 计算两点间距离的平方
 * @param a 第一个点
 * @param b 第二个点
 * @return 返回距离的平方
 */
inline fl vec_distance_sqr(const vec& a, const vec& b) {
    return sqr(a[0] - b[0]) + \
           sqr(a[1] - b[1]) + \
           sqr(a[2] - b[2]);
}

/**
 * @brief 计算向量长度的平方(重载版本)
 * @param v 输入向量
 * @return 返回长度的平方
 */
inline fl sqr(const vec& v) {
    return sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
}

/**
 * @brief 将一个容器的内容追加到另一个容器
 * @tparam T 容器元素类型
 * @param x 目标容器
 * @param y 源容器
 * @return 返回原目标容器的大小
 */
template<typename T>
sz vector_append(std::vector<T>& x, const std::vector<T>& y) {
    sz old_size = x.size();
    x.insert(x.end(), y.begin(), y.end());
    return old_size;
}

/**
 * @brief 查找容器中最小元素的索引
 * @tparam T 容器元素类型
 * @param v 输入容器
 * @return 返回最小元素的索引，空容器返回size()
 */
template<typename T>
sz find_min(const std::vector<T>& v) {
    sz tmp = v.size();
    VINA_FOR_IN(i, v)
        if(i == 0 || v[i] < v[tmp])
            tmp = i;
    return tmp;
}

/**
 * @brief 将角度标准化到[-π, π]范围内
 * @param x 输入角度(引用传递，会被修改)
 */
inline void normalize_angle(fl& x) {
    if     (x >  3*pi) { // 非常大的角度
        fl n = ( x - pi) / (2*pi); // 需要减去多少个2π？
        x -= 2*pi*std::ceil(n); // ceil可能很慢，但这个函数不应该经常调用
        normalize_angle(x);
    }
    else if(x < -3*pi) { // 非常小的角度
        fl n = (-x - pi) / (2*pi);
        x += 2*pi*std::ceil(n);
        normalize_angle(x);
    }
    else if(x >    pi) { // 在(π, 3π]范围内
        x -= 2*pi;
    }
    else if(x <   -pi) { // 在[-3π, -π)范围内
        x += 2*pi;
    }
    assert(x >= -pi && x <= pi);
}

/**
 * @brief 返回标准化后的角度值
 * @param x 输入角度
 * @return 返回[-π, π]范围内的角度值
 */
inline fl normalized_angle(fl x) {
    normalize_angle(x);
    return x;
}

/**
 * @brief 将数值转换为字符串
 * @tparam T 数值类型
 * @param x 输入数值
 * @param width 输出宽度(默认0表示无限制)
 * @param fill 填充字符(默认空格)
 * @return 返回格式化的字符串
 */
template<typename T>
std::string to_string(const T& x, std::streamsize width = 0, char fill = ' ') {
    std::ostringstream out;
    out.fill(fill);
    if(width > 0)
        out << std::setw(width);
    out << x;
    return out.str();
}

/**
 * @brief 计算容器中所有元素的和
 * @tparam T 元素类型
 * @param v 输入容器
 * @return 返回所有元素的和
 */
template<typename T>
T sum(const std::vector<T>& v) {
    T acc = 0;
    VINA_FOR_IN(i, v)
        acc += v[i];
    return acc;
}

/**
 * @brief pK值到自由能的转换因子
 * 
 * 计算公式：
 * K = exp(E/RT)  -- 较低的K和E表示更好的结合
 * pK = -log10(K) => K = 10^(-pK)
 * E = RT ln(K) = RT ln(10^(-pK)) = -RT * ln(10) * pK
 */
const fl pK_to_energy_factor = -8.31 /* RT in J/K/mol */ * 0.001 /* kilo */ * 300 /* K */ / 4.184 /* J/cal */ * std::log(10.0); //  -0.6 kcal/mol * log(10) = -1.38

/**
 * @brief 将pK值转换为自由能(kcal/mol)
 * @param pK 输入pK值
 * @return 返回自由能值
 */
inline fl pK_to_energy(fl pK) { return pK_to_energy_factor * pK; }

/**
 * @brief 打印浮点数
 * @param x 要打印的浮点数
 * @param out 输出流(默认标准输出)
 */
inline void print(fl x, std::ostream& out = std::cout) {
    out << x;
}

/**
 * @brief 打印大小类型数值
 * @param x 要打印的数值
 * @param out 输出流(默认标准输出)
 */
inline void print(sz x, std::ostream& out = std::cout) {
    out << x;
}

/**
 * @brief 打印向量
 * @param v 要打印的向量
 * @param out 输出流(默认标准输出)
 */
inline void print(const vec& v, std::ostream& out = std::cout) {
    out << "(";
    VINA_FOR_IN(i, v) {
        if(i != 0) 
            out << ", ";
        print(v[i], out);
    }
    out << ")";
}

/**
 * @brief 打印容器
 * @tparam T 容器元素类型
 * @param v 要打印的容器
 * @param out 输出流(默认标准输出)
 */
template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
    out << "[";
    VINA_FOR_IN(i, v) {
        if(i != 0) 
            out << " ";
        print(v[i], out);
    }
    out << "]";
}

/**
 * @brief 打印并换行
 * @tparam T 要打印的类型
 * @param x 要打印的对象
 * @param out 输出流(默认标准输出)
 */
template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
    print(x, out);
    out << '\n';
}

/**
 * @brief 检查字符串是否以指定前缀开始
 * @param str 要检查的字符串
 * @param start 前缀字符串
 * @return 如果以指定前缀开始返回true
 */
inline bool starts_with(const std::string& str, const std::string& start) {
    return str.size() >= start.size() && str.substr(0, start.size()) == start;
}

/**
 * @brief 检查容器中是否包含指定元素
 * @tparam T 元素类型
 * @param v 容器
 * @param element 要查找的元素
 * @return 如果包含该元素返回true
 */
template<typename T>
bool has(const std::vector<T>& v, const T& element) {
    return std::find(v.begin(), v.end(), element) != v.end();
}

#endif
