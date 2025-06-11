/**
 * @file triangular_matrix_index.h
 * @brief 三角矩阵索引计算工具
 * 
 * 提供上三角矩阵的一维数组索引计算功能，用于优化对称矩阵的内存存储
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

#ifndef VINA_TRIANGULAR_MATRIX_INDEX_H
#define VINA_TRIANGULAR_MATRIX_INDEX_H

#include "common.h"

/**
 * @brief 计算上三角矩阵在一维数组中的索引位置
 * 
 * 将n×n对称矩阵的上三角部分(包括对角线)压缩存储到一维数组中，
 * 此函数计算矩阵元素(i,j)在压缩数组中的索引位置
 * 
 * @param n 矩阵维度(n×n矩阵)
 * @param i 行索引，必须满足 i <= j
 * @param j 列索引，必须满足 j < n
 * @return 返回元素(i,j)在一维数组中的索引位置
 * 
 * @note 存储顺序为按列优先：
 *       第0列存储元素(0,0)
 *       第1列存储元素(0,1), (1,1)  
 *       第2列存储元素(0,2), (1,2), (2,2)
 *       依此类推...
 *       
 * @note 索引计算公式：index = i + j*(j+1)/2
 *       其中j*(j+1)/2是前j列的元素总数，i是当前列中的偏移量
 */
inline sz triangular_matrix_index(sz n, sz i, sz j) {
    assert(j < n);     // 确保列索引在有效范围内
    assert(i <= j);    // 确保访问的是上三角部分(包括对角线)

    return i + j*(j+1)/2;  // 计算一维数组索引
}

/**
 * @brief 宽松版本的三角矩阵索引计算
 * 
 * 允许i > j的情况，通过交换i和j来确保访问上三角矩阵部分。
 * 这对于对称矩阵特别有用，因为matrix[i][j] = matrix[j][i]
 * 
 * @param n 矩阵维度(n×n矩阵)
 * @param i 行索引
 * @param j 列索引  
 * @return 返回元素(i,j)或(j,i)在一维数组中的索引位置
 * 
 * @note 自动处理索引顺序，确保总是访问上三角部分
 */
inline sz triangular_matrix_index_permissive(sz n, sz i, sz j) {
    return (i <= j) ? triangular_matrix_index(n, i, j)    // 如果i<=j，直接计算
                    : triangular_matrix_index(n, j, i);   // 否则交换i和j后计算
}

#endif
