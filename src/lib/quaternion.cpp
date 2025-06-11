/**
 * @file quaternion.cpp
 * @brief 四元数操作实现文件，用于处理分子旋转变换

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

#include "quaternion.h"

/**
 * @brief 检查四元数是否已归一化
 * @param q 待检查的四元数
 * @return true 如果四元数已归一化
 * @note 用于断言检查，确保四元数模长为1
 */
bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs(q), 1);
}

/**
 * @brief 逐元素比较两个四元数是否近似相等
 * @param a 第一个四元数
 * @param b 第二个四元数
 * @return true 如果所有对应分量都近似相等
 * @note 可能对等价旋转返回false（因为四元数表示旋转的双重性）
 */
bool eq(const qt& a, const qt& b) { // elementwise approximate equality - may return false for equivalent rotations
	return eq(a.R_component_1(), b.R_component_1()) && \
		   eq(a.R_component_2(), b.R_component_2()) && \
		   eq(a.R_component_3(), b.R_component_3()) && \
		   eq(a.R_component_4(), b.R_component_4());
}

/**
 * @brief 从旋转轴和角度创建四元数
 * @param axis 旋转轴（假设为单位向量）
 * @param angle 旋转角度（弧度）
 * @return 对应的四元数
 * @note 使用公式 q = [cos(θ/2), sin(θ/2)*axis]
 */
qt angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
	//assert(eq(tvmet::norm2(axis), 1));
	assert(eq(axis.norm(), 1));
	normalize_angle(angle); // this is probably only necessary if angles can be very big
	fl c = std::cos(angle/2);
	fl s = std::sin(angle/2);
	return qt(c, s*axis[0], s*axis[1], s*axis[2]);
}

/**
 * @brief 从旋转向量创建四元数
 * @param rotation 旋转向量（旋转轴方向 × 旋转角度）
 * @return 对应的四元数
 * @note 旋转向量的模长为旋转角度，方向为旋转轴
 */
qt angle_to_quaternion(const vec& rotation) {
	//fl angle = tvmet::norm2(rotation); 
	fl angle = rotation.norm(); 
	if(angle > epsilon_fl) {
		//vec axis; 
		//axis = rotation / angle;	
		vec axis = (1/angle) * rotation;
		return angle_to_quaternion(axis, angle);
	}
	return qt_identity;
}

/**
 * @brief 将四元数转换为旋转向量表示
 * @param q 输入的归一化四元数
 * @return 旋转向量（轴角表示）
 * @note 旋转向量的模长为旋转角度，方向为旋转轴
 */
vec quaternion_to_angle(const qt& q) {
    assert(quaternion_is_normalized(q));
    const fl c = q.R_component_1(); // 四元数的实部 w = cos(θ/2)
    if(c > -1 && c < 1) { // c理论上应在[-1,1]内，但由于舍入误差可能超出
        fl angle = 2*std::acos(c); // acos返回[0, π]
        if(angle > pi)
            angle -= 2*pi; // 将角度限制在[-π, π]
        vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4()); // 虚部xyz
        fl s = std::sin(angle/2); // 重新计算sin可能不够高效
        if(std::abs(s) < epsilon_fl)
            return zero_vec; // 零旋转情况
        axis *= (angle / s); // 恢复原始旋转向量
        return axis;
    }
    else // 当c = -1或1时，angle/2 = 0或π，因此angle = 0
        return zero_vec;
}

/**
 * @brief 将四元数转换为3x3旋转矩阵
 * @param q 输入的归一化四元数
 * @return 对应的旋转矩阵
 * @note 使用标准的四元数到旋转矩阵转换公式
 */
mat quaternion_to_r3(const qt& q) {
	assert(quaternion_is_normalized(q));

	const fl a = q.R_component_1();
	const fl b = q.R_component_2();
	const fl c = q.R_component_3();
	const fl d = q.R_component_4();

	const fl aa = a*a;
	const fl ab = a*b;
	const fl ac = a*c;
	const fl ad = a*d;
	const fl bb = b*b;
	const fl bc = b*c;
	const fl bd = b*d;
	const fl cc = c*c;
	const fl cd = c*d;
	const fl dd = d*d;

	assert(eq(aa+bb+cc+dd, 1));

	mat tmp;

	// from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
	tmp(0, 0) = (aa+bb-cc-dd);
	tmp(0, 1) = 2*(-ad+bc);
	tmp(0, 2) = 2*(ac+bd);

	tmp(1, 0) = 2*(ad+bc);
	tmp(1, 1) = (aa-bb+cc-dd);
	tmp(1, 2) = 2*(-ab+cd);

	tmp(2, 0) = 2*(-ac+bd);
	tmp(2, 1) = 2*(ab+cd);
	tmp(2, 2) = (aa-bb-cc+dd);

	return tmp;
}

/**
 * @brief 生成随机的单位四元数（均匀分布在旋转空间）
 * @param generator 随机数生成器
 * @return 随机的归一化四元数
 * @note 使用高斯分布生成四个分量，然后归一化
 */
qt random_orientation(rng& generator) {
	qt q(random_normal(0, 1, generator), 
		 random_normal(0, 1, generator), 
		 random_normal(0, 1, generator), 
		 random_normal(0, 1, generator));
	fl nrm = boost::math::abs(q);
	if(nrm > epsilon_fl) {
		q /= nrm;
		assert(quaternion_is_normalized(q));
		return q;
	}
	else 
		return random_orientation(generator); // this call should almost never happen
}

/**
 * @brief 对四元数应用增量旋转
 * @param q 要修改的四元数（引用传递）
 * @param rotation 增量旋转向量
 * @note 计算 q = rotation_quaternion * q，并重新归一化
 */
void quaternion_increment(qt& q, const vec& rotation) {
	assert(quaternion_is_normalized(q));
	q = angle_to_quaternion(rotation) * q;
	quaternion_normalize_approx(q); // normalization added in 1.1.2
	//quaternion_normalize(q); // normalization added in 1.1.2
}

/**
 * @brief 计算从四元数a到四元数b所需的旋转
 * @param b 目标四元数
 * @param a 起始四元数
 * @return 从a转换到b所需的旋转向量
 * @note 计算 b * inv(a) 得到差分旋转
 */
vec quaternion_difference(const qt& b, const qt& a) { // rotation that needs to be applied to convert a to b
	quaternion_is_normalized(a);
	quaternion_is_normalized(b);
	qt tmp = b;
	tmp /= a; // b = tmp * a    =>   b * inv(a) = tmp 
	return quaternion_to_angle(tmp); // already assert normalization
}

/**
 * @brief 将四元数以角度形式打印输出
 * @param q 要打印的四元数
 * @param out 输出流
 */
void print(const qt& q, std::ostream& out) { // print as an angle
	print(quaternion_to_angle(q), out);
}
