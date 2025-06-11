/**
 * @file file.h
 * @brief 文件操作封装类定义
 * @author Dr. Oleg Trott <ot14@columbia.edu>
 * @copyright Copyright (c) 2006-2010, The Scripps Research Institute
 * 
 * 提供基于boost::filesystem的文件流封装类，支持自动错误检测和异常处理
 */
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

#ifndef VINA_FILE_H
#define VINA_FILE_H

#include <boost/filesystem/fstream.hpp>
#include "common.h"

/**
 * @brief 文件操作异常类
 * 
 * 用于表示文件打开或操作失败时的异常信息
 */
struct file_error {
    path name;  ///< 文件路径
    bool in;    ///< 操作类型标志：true表示读取操作，false表示写入操作
    
    /**
     * @brief 构造函数
     * @param name_ 文件路径
     * @param in_ 操作类型标志
     */
    file_error(const path& name_, bool in_) : name(name_), in(in_) {}
};

/**
 * @brief 输入文件流封装类
 * 
 * 继承自boost::filesystem::ifstream，提供自动错误检测功能
 * @warning 不要使用ifstream指针销毁ifile对象，因为没有虚析构函数
 */
struct ifile : public boost::filesystem::ifstream {
    /**
     * @brief 构造函数 - 以默认模式打开文件
     * @param name 文件路径
     * @throws file_error 当文件打开失败时抛出异常
     */
    ifile(const path& name) : boost::filesystem::ifstream(name) {
        if(!(*this))  // 检查文件流状态
            throw file_error(name, true);
    }
    
    /**
     * @brief 构造函数 - 以指定模式打开文件
     * @param name 文件路径
     * @param mode 文件打开模式
     * @throws file_error 当文件打开失败时抛出异常
     */
    ifile(const path& name, std::ios_base::openmode mode) : boost::filesystem::ifstream(name, mode) {
        if(!(*this))  // 检查文件流状态
            throw file_error(name, true);
    }
};

/**
 * @brief 输出文件流封装类
 * 
 * 继承自boost::filesystem::ofstream，提供自动错误检测功能
 * @warning 不要使用ofstream指针销毁ofile对象，因为没有虚析构函数
 */
struct ofile : public boost::filesystem::ofstream {
    /**
     * @brief 构造函数 - 以默认模式打开文件
     * @param name 文件路径
     * @throws file_error 当文件打开失败时抛出异常
     */
    ofile(const path& name) : boost::filesystem::ofstream(name) {
        if(!(*this))  // 检查文件流状态
            throw file_error(name, false);
    }
    
    /**
     * @brief 构造函数 - 以指定模式打开文件
     * @param name 文件路径
     * @param mode 文件打开模式
     * @throws file_error 当文件打开失败时抛出异常
     */
    ofile(const path& name, std::ios_base::openmode mode) : boost::filesystem::ofstream(name, mode) {
        if(!(*this))  // 检查文件流状态
            throw file_error(name, false);
    }
};

#endif
