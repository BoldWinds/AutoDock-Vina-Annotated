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

/// @brief  平滑截断函数
/// @tparam T = fl or vec
/// @param e 能量值
/// @param deriv 导数
/// @param v 截断距离
/// @note 当e<<v时，几乎不变；当e>>v时，能量和导数都显著缩小
template<typename T>
void curl(fl& e, T& deriv, fl v) {
    if(e > 0 && not_max(v)) {
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));  //这个判断是为了数值稳定性
        e *= tmp;
        deriv *= sqr(tmp);
    }
}


inline void curl(fl& e, fl v) {
    if(e > 0 && not_max(v)) {
        fl tmp = (v < epsilon_fl) ? 0 : (v / (v + e));
        e *= tmp;
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
