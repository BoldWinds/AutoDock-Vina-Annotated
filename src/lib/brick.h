/**
 * @file brick.h
 * @brief 空间区域相关的几何计算函数，用于分子对接中的空间距离计算

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

#ifndef VINA_BRICK_H
#define VINA_BRICK_H

#include "common.h"

/**
 * @brief 计算点到线段的最近距离投影点
 * @param begin 线段起始点坐标
 * @param end 线段结束点坐标  
 * @param x 目标点坐标
 * @return fl 返回线段上距离点x最近的坐标值
 * 
 * @note 如果点在线段范围内，返回点本身；否则返回最近的端点
 */
inline fl closest_between(fl begin, fl end, fl x) {
	assert(begin <= end);
	if(x <= begin) return begin;
	else if(x >= end) return end;
	return x;
}

/**
 * @brief 计算三维空间中点到长方体区域的最近点
 * @param begin 长方体区域的最小坐标向量(x_min, y_min, z_min)
 * @param end 长方体区域的最大坐标向量(x_max, y_max, z_max)
 * @param v 目标点的坐标向量
 * @return vec 返回长方体表面或内部距离目标点最近的点坐标
 * 
 * @note 对每个坐标轴分别应用closest_between函数，确保结果点在长方体范围内
 */
inline vec brick_closest(const vec& begin, const vec& end, const vec& v) {
	vec tmp;
	VINA_FOR_IN(i, tmp)
		tmp[i] = closest_between(begin[i], end[i], v[i]);
	return tmp;
}

/**
 * @brief 计算点到长方体区域的最短距离的平方
 * @param begin 长方体区域的最小坐标向量
 * @param end 长方体区域的最大坐标向量
 * @param v 目标点的坐标向量
 * @return fl 返回点到长方体的最短距离的平方值
 * 
 * @note 先找到长方体上最近的点，然后计算该点与目标点的距离平方
 * @note 返回距离平方是为了避免开平方运算，提高计算效率
 */
inline fl brick_distance_sqr(const vec& begin, const vec& end, const vec& v) {
	vec closest; closest = brick_closest(begin, end, v);
	return vec_distance_sqr(closest, v);
}

#endif
