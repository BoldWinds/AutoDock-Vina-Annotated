/**
 * @file cache.h
 * @brief 管理所有原子类型的能量网格，实现能量计算接口

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

#ifndef VINA_CACHE_H
#define VINA_CACHE_H

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/static_assert.hpp>
#include "igrid.h"
#include "grid.h"
#include "model.h"
#include "file.h"
#include "szv_grid.h"


struct precalculate;
struct model;

/**
 * @brief 管理所有原子类型的能量网格，提供统一的能量计算接口
 */
struct cache : public igrid {
private:
	grid_dims m_gd;
	fl m_slope; // does not get (de-)serialized
	std::vector<grid> m_grids;      // 存储所有原子类型的网格
public:
    cache(fl slope=1e6) : m_slope(slope), m_grids(XS_TYPE_SIZE) {}
	cache(const grid_dims& gd, fl slope=1e6) : m_gd(gd), m_slope(slope), m_grids(XS_TYPE_SIZE) {}
	fl eval      (const model& m, fl v) const; // needs m.coords // clean up
	fl eval_intra(      model& m, fl v) const; // needs m.coords // clean up
	fl eval_deriv(      model& m, fl v) const; // needs m.coords, sets m.minus_forces // clean up
    grid_dims gd() const { return m_gd; }
    vec corner1() const { vec corner(m_gd[0].begin, m_gd[1].begin, m_gd[2].begin); return corner; }
    vec corner2() const { vec corner(m_gd[0].end, m_gd[1].end, m_gd[2].end); return corner; }
    bool is_in_grid(const model &m, fl margin=0.0001) const;
    bool is_atom_type_grid_initialized(sz t) const { return m_grids[t].initialized(); }
    bool are_atom_types_grid_initialized(szv atom_types) const;
    void read(const std::string& str);
    void write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename="NULL",
               const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");
    /**
    * @brief 填充网格数据(m_grids)
    * @param m 模型对象
    * @param p 原子类型的预计算对象
    * @param atom_types_needed 需要填充的原子类型列表
    */
	void populate(const model& m, const precalculate& p, const szv& atom_types_needed);
};

#endif
