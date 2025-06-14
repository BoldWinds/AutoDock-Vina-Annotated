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

#ifndef VINA_H
#define VINA_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include <boost/algorithm/string.hpp>
//#include <openbabel/mol.h>
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "conf.h"
#include "model.h"
#include "common.h"
#include "cache.h"
#include "non_cache.h"
#include "ad4cache.h"
#include "quasi_newton.h"
#include "coords.h" // add_to_output_container
#include "utils.h"
#include "scoring_function.h"
#include "precalculate.h"
#include <memory>

/**
 * @class vina_runtime_error
 * @brief Vina运行时异常类
 * 
 * 继承自std::exception，用于处理Vina运行过程中的错误情况
 */
class vina_runtime_error : public std::exception {
public:
    /**
     * @brief 构造函数
     * @param message 错误信息字符串
     */
    explicit vina_runtime_error(const std::string & message)
            : m_message("\n\nVina runtime error: " + message + "\n") {}

    /**
     * @brief 获取异常信息
     * @return 异常信息的C风格字符串
     */
    virtual const char* what() const throw () {
        return m_message.c_str();
    }

private:
    const std::string m_message;
};

/**
 * @class Vina
 * @brief AutoDock Vina主要接口类
 * 
 * 该类封装了AutoDock Vina分子对接的全部功能，包括：
 * - 受体和配体的设置
 * - 评分函数的配置
 * - 网格地图的计算
 * - 分子对接的全局搜索
 * - 结果的优化和输出
 */
class Vina {
private:
    // 模型和构象相关成员变量
    model m_receptor;                   ///< 受体分子模型
    model m_model;                      ///< 配体分子模型
    output_container m_poses;           ///< 对接结果构象容器
    //OpenBabel::OBMol m_mol;
    bool m_receptor_initialized;        ///< 受体是否已初始化标志
    bool m_ligand_initialized;          ///< 配体是否已初始化标志

    // 评分函数相关成员变量
    scoring_function_choice m_sf_choice;        ///< 评分函数选择
    flv m_weights;                             ///< 评分函数权重向量
    std::shared_ptr<ScoringFunction> m_scoring_function;    ///< 评分函数对象指针
    precalculate_byatom m_precalculated_byatom;            ///< 按原子类型预计算数据
    precalculate m_precalculated_sf;                       ///< 评分函数预计算数据

    // 网格地图相关成员变量
    cache m_grid;                       ///< Vina网格缓存
    ad4cache m_ad4grid;                 ///< AutoDock4网格缓存
    non_cache m_non_cache;              ///< 非缓存计算对象
    bool m_map_initialized;             ///< 网格地图是否已初始化标志

    // 全局搜索相关成员变量
    int m_cpu;                          ///< CPU核心数
    int m_seed;                         ///< 随机数种子

    // 其他成员变量
    int m_verbosity;                                    ///< 详细程度等级
    bool m_no_refine;                                   ///< 是否禁用精化标志
    std::function<void(double)>* m_progress_callback;   ///< 进度回调函数指针

    /**
     * @brief 生成Vina格式的备注信息
     * @param pose 构象输出类型引用
     * @param lb 下界值
     * @param ub 上界值
     * @return 格式化的备注字符串
     */
	std::string vina_remarks(output_type& pose, fl lb, fl ub);

    /**
     * @brief 移除冗余构象
     * @param in 输入构象容器
     * @param min_rmsd 最小RMSD阈值
     * @return 去除冗余后的构象容器
     * 
     * @note 通过RMSD比较移除相似的构象，保留多样性
     */
	output_container remove_redundant(const output_container& in, fl min_rmsd);

    /**
     * @brief 设置力场参数
     */
	void set_forcefield();

    /**
     * @brief 优化指定构象（私有重载版本）
     * @param out 要优化的构象
     * @param max_steps 最大优化步数（0表示使用默认值）
     * @return 优化后的能量信息向量
     */
	std::vector<double> optimize(output_type& out, const int max_steps=0);

    /**
     * @brief 生成随机数种子
     * @param seed 输入种子（0表示使用当前时间）
     * @return 生成的随机数种子
     */
	int generate_seed(const int seed=0);

public:
    /**
     * @brief 构造函数
     * @param sf_name 评分函数名称（默认为"vina"，可选"vinardo"或"ad4"）
     * @param cpu CPU核心数（0表示自动检测）
     * @param seed 随机数种子（0表示使用当前时间）
     * @param verbosity 详细程度等级（默认为1）
     * @param no_refine 是否禁用精化步骤（默认为false）
     * @param progress_callback 进度回调函数指针（默认为NULL）
     * 
     * @note 构造函数会根据sf_name选择相应的评分函数并设置默认权重
     */
	Vina(const std::string &sf_name="vina", int cpu=0, int seed=0, int verbosity=1, bool no_refine=false, std::function<void(double)>* progress_callback = NULL) {
		m_verbosity = verbosity;
		m_receptor_initialized = false;
		m_ligand_initialized = false;
		m_map_initialized = false;
		m_seed = generate_seed(seed);
		m_no_refine = no_refine;
		m_progress_callback = progress_callback;

		// 查找CPU核心数
		if (cpu <= 0) {
			unsigned num_cpus = boost::thread::hardware_concurrency();

			if (num_cpus > 0) {
				m_cpu = num_cpus;
			} else {
				std::cerr << "WARNING: Could not determined the number of concurrent thread supported on this machine. ";
				std::cerr << "You might need to set it manually using cpu argument or fix the issue.\n";
				exit(EXIT_FAILURE);
			}
		} else {
			m_cpu = cpu;
		}
		
		// 根据评分函数名称设置相应的评分函数和权重
		if (sf_name.compare("vina") == 0) {
			m_sf_choice = SF_VINA;
			set_vina_weights();
		} else if (sf_name.compare("vinardo") == 0) {
			m_sf_choice = SF_VINARDO;
			set_vinardo_weights();
		} else if (sf_name.compare("ad4") == 0) {
			m_sf_choice = SF_AD42;
			set_ad4_weights();
		} else {
			std::cerr << "ERROR: Scoring function " << sf_name << " not implemented (choices: vina, vinardo or ad4)\n";
			exit (EXIT_FAILURE);
		}
	}
	// Destructor
	~Vina();

    /**
     * @brief 显示引用信息
     */
	void cite();

    /**
     * @brief 获取随机数种子
     * @return 当前使用的随机数种子
     */
	int seed() { return m_seed; }

    /**
     * @brief 设置受体分子
     * @param rigid_name 刚性受体文件名（PDBQT格式）
     * @param flex_name 柔性受体文件名（PDBQT格式，可选）
     */
	void set_receptor(const std::string &rigid_name=std::string(), const std::string &flex_name=std::string());

    /**
     * @brief 从字符串设置配体分子
     * @param ligand_string 包含配体信息的PDBQT格式字符串
     */
	void set_ligand_from_string(const std::string &ligand_string);

    /**
     * @brief 从字符串向量设置多个配体分子
     * @param ligand_string 包含多个配体信息的PDBQT格式字符串向量
     */
	void set_ligand_from_string(const std::vector<std::string> &ligand_string);

    /**
     * @brief 从文件设置配体分子
     * @param ligand_name 配体文件名（PDBQT格式）
     */
	void set_ligand_from_file(const std::string& ligand_name);

    /**
     * @brief 从文件向量设置多个配体分子
     * @param ligand_name 配体文件名向量（PDBQT格式）
     */
	void set_ligand_from_file(const std::vector<std::string>& ligand_name);

	//void set_ligand(OpenBabel::OBMol* mol);
	//void set_ligand(std::vector<OpenBabel::OBMol*> mol);

    /**
     * @brief 设置Vina评分函数的权重参数
     * @param weight_gauss1 高斯项1权重（默认-0.035579）
     * @param weight_gauss2 高斯项2权重（默认-0.005156）
     * @param weight_repulsion 排斥项权重（默认0.840245）
     * @param weight_hydrophobic 疏水项权重（默认-0.035069）
     * @param weight_hydrogen 氢键项权重（默认-0.587439）
     * @param weight_glue 胶水项权重（默认50）
     * @param weight_rot 旋转项权重（默认0.05846）
     */
	void set_vina_weights(double weight_gauss1=-0.035579, double weight_gauss2=-0.005156,
						       double weight_repulsion=0.840245, double weight_hydrophobic=-0.035069,
						       double weight_hydrogen=-0.587439, double weight_glue=50,
						       double weight_rot=0.05846);

    /**
     * @brief 设置Vinardo评分函数的权重参数
     * @param weight_gauss1 高斯项权重（默认-0.045）
     * @param weight_repulsion 排斥项权重（默认0.8）
     * @param weight_hydrophobic 疏水项权重（默认-0.035）
     * @param weight_hydrogen 氢键项权重（默认-0.600）
     * @param weight_glue 胶水项权重（默认50）
     * @param weight_rot 旋转项权重（默认0.05846）
     */
	void set_vinardo_weights(double weight_gauss1=-0.045,
							       double weight_repulsion=0.8, double weight_hydrophobic=-0.035,
							       double weight_hydrogen=-0.600, double weight_glue=50,
							       double weight_rot=0.05846);

    /**
     * @brief 设置AutoDock4评分函数的权重参数
     * @param weight_ad4_vdw 范德华项权重（默认0.1662）
     * @param weight_ad4_hb 氢键项权重（默认0.1209）
     * @param weight_ad4_elec 静电项权重（默认0.1406）
     * @param weight_ad4_dsolv 去溶剂化项权重（默认0.1322）
     * @param weight_glue 胶水项权重（默认50）
     * @param weight_ad4_rot 旋转项权重（默认0.2983）
     */
	void set_ad4_weights(double weight_ad4_vdw=0.1662, double weight_ad4_hb=0.1209,
						      double weight_ad4_elec=0.1406, double weight_ad4_dsolv=0.1322,
						      double weight_glue=50, double weight_ad4_rot=0.2983);

    /**
     * @brief 根据配体分子计算网格维度
     * @param buffer_size 缓冲区大小（默认4Å）
     * @return 包含网格尺寸信息的向量
     * 
     * @note 自动计算包围配体的最小网格盒子，并添加指定的缓冲区
     */
	std::vector<double> grid_dimensions_from_ligand(double buffer_size=4);

    /**
     * @brief 计算Vina评分函数的网格地图
     * @param center_x 网格中心X坐标
     * @param center_y 网格中心Y坐标
     * @param center_z 网格中心Z坐标
     * @param size_x 网格X方向尺寸
     * @param size_y 网格Y方向尺寸
     * @param size_z 网格Z方向尺寸
     * @param granularity 网格粒度（默认0.5Å）
     * @param force_even_voxels 是否强制偶数体素（默认false）
     * 
     * @note 预计算评分函数在网格点上的值，用于加速对接计算
     */
	void compute_vina_maps(double center_x, double center_y, double center_z,
								  double size_x, double size_y, double size_z,
								  double granularity=0.5, bool force_even_voxels=false);

    /**
     * @brief 从文件加载预计算的网格地图
     * @param maps 网格地图文件路径
     */
	void load_maps(std::string maps);

    /**
     * @brief 随机化配体构象
     * @param max_steps 最大随机化步数（默认10000）
     * 
     * @note 生成随机的配体起始构象，用于全局搜索
     */
	void randomize(const int max_steps=10000);

    /**
     * @brief 计算当前配体构象的评分
     * @return 包含总能量和各项能量贡献的向量
     */
	std::vector<double> score();

    /**
     * @brief 计算当前配体构象的评分（指定分子内能量）
     * @param intramolecular_energy 分子内能量
     * @return 包含总能量和各项能量贡献的向量
     */
	std::vector<double> score(double intramolecular_energy);

    /**
     * @brief 优化当前配体构象
     * @param max_steps 最大优化步数（0表示使用默认值）
     * @return 优化后的能量信息向量
     * 
     * @note 使用准牛顿法进行局部优化
     */
	std::vector<double> optimize(const int max_steps=0);

    /**
     * @brief 执行全局搜索算法
     * @param exhaustiveness 搜索详尽程度（默认8）
     * @param n_poses 保留的构象数量（默认20）
     * @param min_rmsd 最小RMSD阈值（默认1.0Å）
     * @param max_evals 最大评估次数（0表示无限制）
     * 
     * @note 使用蒙特卡洛搜索算法寻找最优结合构象
     */
	void global_search(const int exhaustiveness=8, const int n_poses=20, const double min_rmsd=1.0, const int max_evals=0);

    /**
     * @brief 获取构象结果的PDBQT格式字符串
     * @param how_many 返回的构象数量（默认9）
     * @param energy_range 能量范围阈值（默认3.0 kcal/mol）
     * @return 包含构象信息的PDBQT格式字符串
     */
	std::string get_poses(int how_many=9, double energy_range=3.0);

    /**
     * @brief 获取构象的坐标信息
     * @param how_many 返回的构象数量（默认9）
     * @param energy_range 能量范围阈值（默认3.0 kcal/mol）
     * @return 包含各构象坐标的二维向量
     */
	std::vector< std::vector<double> > get_poses_coordinates(int how_many=9, double energy_range=3.0);

    /**
     * @brief 获取构象的能量信息
     * @param how_many 返回的构象数量（默认9）
     * @param energy_range 能量范围阈值（默认3.0 kcal/mol）
     * @return 包含各构象能量的二维向量
     */
	std::vector< std::vector<double> > get_poses_energies(int how_many=9, double energy_range=3.0);

    /**
     * @brief 写入单个构象到文件
     * @param output_name 输出文件名
     * @param remark 备注信息（可选）
     */
	void write_pose(const std::string &output_name, const std::string &remark = std::string());

    /**
     * @brief 写入多个构象到文件
     * @param output_name 输出文件名
     * @param how_many 写入的构象数量（默认9）
     * @param energy_range 能量范围阈值（默认3.0 kcal/mol）
     */
	void write_poses(const std::string &output_name, int how_many=9, double energy_range=3.0);

    /**
     * @brief 写入网格地图文件
     * @param map_prefix 地图文件前缀（默认"receptor"）
     * @param gpf_filename GPF文件名（默认"NULL"）
     * @param fld_filename FLD文件名（默认"NULL"）
     * @param receptor_filename 受体文件名（默认"NULL"）
     */
	void write_maps(const std::string& map_prefix="receptor", const std::string& gpf_filename="NULL",
					    const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");

    /**
     * @brief 显示评分信息
     * @param energies 能量向量
     */
	void show_score(const std::vector<double> energies);
};

#endif
