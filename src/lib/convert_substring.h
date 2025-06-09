/**
 * @file convert_substring.h
 * @brief 字符串子串转换工具函数，用于AutoDock Vina中的数据解析
 *
 * Copyright (c) 2006-2010, The Scripps Research Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 * Author: Dr. Oleg Trott <ot14@columbia.edu>, 
 *         The Olson Lab, 
 *         The Scripps Research Institute
 */

#ifndef VINA_CONVERT_SUBSTRING_H
#define VINA_CONVERT_SUBSTRING_H

#include <cctype> // for isspace
#include <boost/lexical_cast.hpp>
#include "common.h"

/**
 * @brief 转换失败时抛出的异常结构体
 */
struct bad_conversion {};

/**
 * @brief 将字符串的指定子串转换为目标类型T
 * @tparam T 目标转换类型
 * @param str 源字符串
 * @param i 起始索引（1基索引）
 * @param j 结束索引（1基索引，包含）
 * @return 转换后的T类型值
 * @throws bad_conversion 当索引无效或转换失败时抛出
 * @note 使用1基索引系统，会自动跳过前导空白字符
 */
template<typename T>
T convert_substring(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
    // 检查索引边界：i必须>=1，i不能超过j+1，j不能超过字符串长度
    if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();

    // 跳过前导空白字符
    while(i <= j && std::isspace(str[i-1]))
        ++i;

    T tmp;
    try {
        // 使用boost::lexical_cast进行类型转换
        // 转换为0基索引：substr(i-1, j-i+1)
        tmp = boost::lexical_cast<T>(str.substr(i-1, j-i+1));
    }
    catch(...) {
        // 捕获所有boost::lexical_cast可能抛出的异常
        throw bad_conversion();
    }
    return tmp;
}

/**
 * @brief 检查指定子串是否只包含空白字符
 * @param str 源字符串
 * @param i 起始索引（1基索引）
 * @param j 结束索引（1基索引，包含）
 * @return true 如果子串只包含空白字符，false 否则
 * @throws bad_conversion 当索引无效时抛出
 */
inline bool substring_is_blank(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
    // 检查索引边界
    if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();
    // 遍历指定范围内的所有字符
    VINA_RANGE(k, i-1, j)
        if(!std::isspace(str[k]))  // 发现非空白字符立即返回false
            return false;
    return true;  // 所有字符都是空白字符
}

/**
 * @brief convert_substring针对unsigned类型的特化版本
 * @param str 源字符串
 * @param i 起始索引（1基索引）
 * @param j 结束索引（1基索引，包含）
 * @return 转换后的unsigned值
 * @throws bad_conversion 当索引无效、转换失败或结果为负数时抛出
 * @note 先转换为int类型，然后检查是否为负数，这是为了处理早期boost::lexical_cast
 *       对unsigned类型处理"-123"等负数字符串的问题
 */
template<>
inline unsigned convert_substring<unsigned>(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
    int tmp = convert_substring<int>(str, i, j);  // 先转换为int类型
    if(tmp < 0) throw bad_conversion();  // 检查是否为负数
    return static_cast<unsigned>(tmp);  // 安全转换为unsigned
}

#endif
