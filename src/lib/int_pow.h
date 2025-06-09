/**
 * @file int_pow.h
 * @brief 编译时整数幂运算模板函数实现
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

#ifndef VINA_INT_POW_H
#define VINA_INT_POW_H

#include "common.h" // for fl

/**
 * @brief 编译时整数幂运算模板函数
 * 
 * 使用递归模板技术计算 x 的 n 次幂。编译器会在编译时将递归完全展开，
 * 生成一系列直接乘法运算，避免运行时循环开销。
 * 
 * @tparam n 指数（必须是编译时常量）
 * @param x 底数（浮点数类型）
 * @return fl 返回 x^n 的计算结果
 * 
 * @note 递归展开示例：int_pow<3>(x) → ((1*x)*x)*x
 * @note 编译器优化良好，即使对于较大的 n 值也能生成高效代码
 */
template<unsigned n>
inline fl int_pow(fl x) {
    return int_pow<n-1>(x)*x;
}

/**
 * @brief 模板特化：零次幂的情况
 * 
 * 任何数的零次幂都等于1，这是递归模板的终止条件。
 * 
 * @param x 底数（未使用，任何数的0次幂都是1）
 * @return fl 返回1
 * 
 * @note 数学基础：x^0 = 1（x ≠ 0时）
 */
template<>
inline fl int_pow<0>(fl x) {
    return 1; // 递归终止条件：任何数的0次幂等于1
}

#endif
