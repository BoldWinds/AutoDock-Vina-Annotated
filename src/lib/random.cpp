/**
 * @file random.cpp
 * @brief 随机数生成相关函数实现

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

#include <ctime> // for time (for seeding)

#include "random.h"

/**
 * @brief 生成指定范围内的随机浮点数
 * @param a 范围下界（不包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回[a, b]范围内的随机浮点数
 * 
 * 使用Boost库的uniform_real分布生成均匀分布的随机浮点数。
 * 函数会进行断言检查以确保参数有效性和返回值的正确性。
 */
fl random_fl(fl a, fl b, rng& generator) { // expects a < b, returns rand in [a, b]
    assert(a < b); // BOOST also asserts a < b
    typedef boost::uniform_real<fl> distr;
    boost::variate_generator<rng&, distr> r(generator, distr(a, b));
    fl tmp = r();
    assert(tmp >= a);
    assert(tmp <= b);
    return tmp;
}

/**
 * @brief 生成符合正态分布的随机浮点数
 * @param mean 正态分布的均值
 * @param sigma 正态分布的标准差
 * @param generator 随机数生成器引用
 * @return 返回符合指定均值和标准差的正态分布随机浮点数
 * 
 * 使用Boost库的normal_distribution生成正态分布（高斯分布）的随机数。
 */
fl random_normal(fl mean, fl sigma, rng& generator) { // expects sigma >= 0
    assert(sigma >= 0); // BOOST asserts this as well
    typedef boost::normal_distribution<fl> distr;
    boost::variate_generator<rng&, distr> r(generator, distr(mean, sigma));
    return r();
}

/**
 * @brief 生成指定范围内的随机整数
 * @param a 范围下界（包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回[a, b]范围内的随机整数
 * 
 * 使用Boost库的uniform_int分布生成均匀分布的随机整数。
 * 函数会进行断言检查以确保参数有效性和返回值的正确性。
 */
int random_int(int a, int b, rng& generator) { // expects a <= b, returns rand in [a, b]
    assert(a <= b); // BOOST asserts this as well
    typedef boost::uniform_int<int> distr;
    boost::variate_generator<rng&, distr> r(generator, distr(a, b));
    int tmp = r();
    assert(tmp >= a);
    assert(tmp <= b);
    return tmp;
}

/**
 * @brief 生成指定范围内的随机无符号整数
 * @param a 范围下界（包含）
 * @param b 范围上界（包含）
 * @param generator 随机数生成器引用
 * @return 返回[a, b]范围内的随机无符号整数
 * 
 * 通过调用random_int函数生成整数，然后转换为无符号类型。
 * 包含多重断言检查以确保类型转换的安全性。
 */
sz random_sz(sz a, sz b, rng& generator) { // expects a <= b, returns rand in [a, b]
    assert(a <= b);
    assert(int(a) >= 0);
    assert(int(b) >= 0);
    int i = random_int(int(a), int(b), generator);
    assert(i >= 0);
    assert(i >= int(a));
    assert(i <= int(b));
    return static_cast<sz>(i);
}

/**
 * @brief 生成单位球体内的随机三维向量
 * @param generator 随机数生成器引用
 * @return 返回以原点为中心、半径为1的球体内的随机三维向量
 * 
 * 使用拒绝采样方法：在[-1,1]³的立方体内随机生成点，
 * 只接受距离原点小于1的点。平均需要运行约2次。
 */
vec random_inside_sphere(rng& generator) {
    while(true) { // on average, this will have to be run about twice
        fl r1 = random_fl(-1, 1, generator);
        fl r2 = random_fl(-1, 1, generator);
        fl r3 = random_fl(-1, 1, generator);

        vec tmp(r1, r2, r3);
        if(sqr(tmp) < 1)  // 检查向量的平方长度是否小于1
            return tmp;
    }
}

/**
 * @brief 生成指定矩形框内的随机三维向量
 * @param corner1 矩形框的第一个角点坐标
 * @param corner2 矩形框的第二个角点坐标
 * @param generator 随机数生成器引用
 * @return 返回矩形框内的随机三维向量
 * 
 * 在每个维度上独立生成[corner1[i], corner2[i]]范围内的随机坐标。
 */
vec random_in_box(const vec& corner1, const vec& corner2, rng& generator) { // expects corner1[i] < corner2[i]
    vec tmp;
    VINA_FOR_IN(i, tmp)  // 对向量的每个维度进行迭代
        tmp[i] = random_fl(corner1[i], corner2[i], generator);
    return tmp;
}

/**
 * @brief 自动生成随机种子
 * @return 返回基于系统熵的随机种子值
 * 
 * 使用现代C++的随机设备获取高质量的随机种子：
 * 1. 使用std::random_device从操作系统熵池获取随机数
 * 2. 用该随机数初始化64位梅森旋转算法生成器
 * 3. 通过均匀分布生成最终的种子值
 * 
 * 这种方法比基于PID和时间的传统方法更加安全可靠。
 */
int auto_seed() {
    // Seed generator, fix previous seed generator based on PID and time
    // Source: https://stackoverflow.com/questions/22883840/c-get-random-number-from-0-to-max-long-long-integer
    std::random_device rd; // Get a random seed from the OS entropy device, or whatever
    std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937 generator
                               // and seed it with entropy.

    // Define the distribution, by default it goes from 0 to MAX(unsigned int)
    // or what have you.
    std::uniform_int_distribution<unsigned int> distr;
    return distr(eng);
}
