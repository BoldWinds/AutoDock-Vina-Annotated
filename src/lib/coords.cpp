/**
 * @file coords.cpp
 * @brief 实现计算RMSD上界计算、查找最相似构象和添加对接结果到输出容器

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

#include "coords.h"

/**
 * @brief 计算两组坐标向量间的RMSD上界
 * @param a 第一组坐标向量
 * @param b 第二组坐标向量
 * @return 计算得到的RMSD值
 * @note 使用平方根均值方差公式计算分子构象间的几何差异
 */
fl rmsd_upper_bound(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a) 
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}

/**
 * @brief 在输出容器中查找与给定坐标最相似的构象
 * @param a 待比较的坐标向量
 * @param b 输出容器，包含多个构象
 * @return pair<索引, RMSD值> 最相似构象的索引和对应的RMSD距离
 */
std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
	std::pair<sz, fl> tmp(b.size(), max_fl);
	VINA_FOR_IN(i, b) {
		fl res = rmsd_upper_bound(a, b[i].coords);
		if(i == 0 || res < tmp.second)
			tmp = std::pair<sz, fl>(i, res);
	}
	return tmp;
}

/**
 * @brief 向输出容器添加新的对接结果，并检查rmsd是否符合要求、去重，并按照结合能排序
 */
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size) {
	std::pair<sz, fl> closest_rmsd = find_closest(t.coords, out);
	if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
		if(t.e < out[closest_rmsd.first].e) { // the new one is better, apparently
			out[closest_rmsd.first] = t; // FIXME? slow
		}
	}
	else { // nothing similar
		if(out.size() < max_size)
			out.push_back(new output_type(t)); // the last one had the worst energy - replacing 
		else
			if(!out.empty() && t.e < out.back().e) // FIXME? - just changed
				out.back() = t; // FIXME? slow
	}
	out.sort();
}
