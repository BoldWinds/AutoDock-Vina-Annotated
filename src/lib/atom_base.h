/**
 * @file atom_base.h
 * @brief 给原子添加电荷

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

#ifndef VINA_ATOM_BASE_H
#define VINA_ATOM_BASE_H

#include "atom_type.h"

/**
 * @brief 原子基础类
 */
struct atom_base : public atom_type {
    fl charge;  ///< 原子电荷(原子单位)

    /**
     * @brief 默认构造函数
     * 
     * 初始化电荷为0，原子类型通过父类构造函数初始化为未分配状态
     */
    atom_base() : charge(0) {}

private:
    friend class boost::serialization::access;

    /**
     * @brief 序列化函数，用于boost序列化支持
     * @param ar 序列化归档对象
     * @param version 版本号(未使用)
     * @note 先序列化父类atom_type，再序列化自身的charge属性
     */
    template<class Archive> 
    void serialize(Archive& ar, const unsigned version) {
        ar & boost::serialization::base_object<atom_type>(*this);
        ar & charge;
    }
};

#endif
