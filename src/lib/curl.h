/**
 * @file curl.h
 * @brief 能量截断函数实现，用于AutoDock Vina分子对接中的能量平滑处理

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

#ifndef VINA_CURL_H
#define VINA_CURL_H

#include "common.h"

#if 1 // 优先使用软截断而非硬截断
/**
 * @brief 软截断函数模板 - 对能量值进行平滑截断并计算导数
 * @tparam T 导数类型，可以是fl（标量）或vec（向量）
 * @param[in,out] e 输入/输出能量值，会被修改为截断后的值
 * @param[in,out] deriv 输入/输出导数值，会被修改为截断后的导数
 * @param[in] v 截断阈值参数
 * 
 * @note 使用平滑函数 e' = e * (v/(v+e))，保证梯度连续性
 * @note 当v很小时（< epsilon_fl），tmp设为0以避免数值不稳定
 * @note 导数按链式法则计算：deriv' = deriv * (v/(v+e))^2
 */
template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
    if(e > 0 && not_max(v)) { // 仅当能量为正且v不是最大值时进行截断
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e)); // 计算截断因子，避免除零
        e *= tmp;          // 应用截断到能量值
        deriv *= sqr(tmp); // 应用截断到导数（平方项来自链式法则）
    }
}

/**
 * @brief 软截断函数 - 仅对能量值进行平滑截断，不涉及导数计算
 * @param[in,out] e 输入/输出能量值，会被修改为截断后的值
 * @param[in] v 截断阈值参数
 * 
 * @note 与模板版本相同的平滑函数，但不计算导数
 */
inline void curl(fl& e, fl v) {
    if(e > 0 && not_max(v)) { // 仅当能量为正且v不是最大值时进行截断
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e)); // 计算截断因子
        e *= tmp; // 应用截断到能量值
    }
}

#else

/**
 * @brief 硬截断函数模板 - 对能量值进行直接阈值截断
 * @tparam T 导数类型，可以是fl（标量）或vec（向量）
 * @param[in,out] e 输入/输出能量值，如果超过阈值则设为阈值
 * @param[in,out] deriv 输入/输出导数值，截断时设为0
 * @param[in] v 截断阈值
 * 
 * @note 硬截断：如果e > v，则直接设e = v，导数设为0
 * @note 这种方法简单但在截断点处梯度不连续
 */
template<typename T> // T = fl or vec
void curl(fl& e, T& deriv, fl v) {
    if(e > v) {
        e = v;     // 直接截断到阈值
        deriv = 0; // 导数置零
    }
}

/**
 * @brief 硬截断函数 - 仅对能量值进行直接阈值截断
 * @param[in,out] e 输入/输出能量值，如果超过阈值则设为阈值
 * @param[in] v 截断阈值
 */
inline void curl(fl& e, fl v) {
    if(e > v) {
        e = v; // 直接截断到阈值
    }
}
#endif

#endif
