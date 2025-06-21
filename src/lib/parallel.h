/**
 * @file parallel.h
 * @brief AutoDock Vina多线程并行计算框架头文件

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

#ifndef VINA_PARALLEL_H
#define VINA_PARALLEL_H

#include <vector>

#include "common.h"

#include <boost/optional.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp> 


/**
 * @brief 并行for循环模板类
 * @tparam F 函数对象类型，用于处理每个索引
 * @tparam Sync 同步模式标志，false为异步模式，true为同步模式，只用到了同步模式
 */
template<typename F, bool Sync = false>
struct parallel_for : private boost::thread_group {
    // 构造函数，创建指定数量的工作线程
	parallel_for(const F* f, sz num_threads) : m_f(f), destructing(false), size(0), thread_finished(num_threads, true), count_finished(0), num_threads(num_threads) {
        VINA_FOR(i, num_threads)
            create_thread(aux(i, this));
    }
    /**
     * @brief 执行并行计算，阻塞直到所有线程完成处理
     */
	void run(sz size_) {
		boost::mutex::scoped_lock self_lk(self);
        size = size_;
        count_finished = 0;
		VINA_FOR_IN(i, thread_finished) 
			thread_finished[i] = false;
        cond.notify_all(); // many things modified
        while(count_finished < num_threads) // wait until processing of all elements is thread_finished
            busy.wait(self_lk);
    }
    /**
     * @brief 析构函数，安全关闭所有工作线程
     */
    virtual ~parallel_for() {
        {
			boost::mutex::scoped_lock self_lk(self);
            destructing = true;
            cond.notify_all(); // destructing modified
        }
        join_all(); 
    }
private:
    /**
     * @brief 工作线程的主循环函数
     * @param offset 线程偏移量，用于确定处理的起始位置
     */
	void loop(sz offset) {
		while(boost::optional<sz> sz_option = get_size(offset)) {
			sz s = sz_option.get();
			for(sz i = offset; i < s; i += num_threads)
				(*m_f)(i);
			{
				boost::mutex::scoped_lock self_lk(self);
				thread_finished[offset] = true;
				++count_finished;
				busy.notify_one();
			}
		}
	}
    /**
     * @brief 线程辅助结构体，用于启动工作线程
     */
    struct aux {
		sz offset;
        parallel_for* par;
		aux(sz offset, parallel_for* par) : offset(offset), par(par) {}
		void operator()() const { par->loop(offset); }
    };
    const F* m_f;                    ///< 函数对象指针（不持有本地副本）
    boost::condition cond;           ///< 条件变量，用于同步工作线程
    boost::condition busy;           ///< 条件变量，用于通知主线程
    bool destructing;                ///< 析构标志
    sz size;                         ///< 传递给run()的向量大小
    std::vector<bool> thread_finished; ///< 线程完成状态数组
    sz count_finished;               ///< 已完成的线程数量
    sz num_threads;                  ///< 工作线程数量，也是步长
    boost::mutex self;               ///< 互斥锁，保护可变成员的读写
    
    /**
     * @brief 获取当前任务大小（线程安全）
     * @param offset 线程偏移量
     * @return 任务大小的可选值，如果正在析构则返回空
     */
	boost::optional<sz> get_size(sz offset) {
		boost::mutex::scoped_lock self_lk(self);
        while(!destructing && thread_finished[offset])
            cond.wait(self_lk);
		if(destructing) return boost::optional<sz>(); // wrap it up!
        return size;
    }
};

/**
 * @brief 并行for循环模板类的同步模式特化版本，也是在实际代码中唯一使用的版本
 * @tparam F 函数对象类型
 */
template<typename F>
struct parallel_for<F, true> : private boost::thread_group {
	parallel_for(const F* f, sz num_threads) : m_f(f), destructing(false), size(0), started(0), finished(0) {
		a.par = this; // VC8 warning workaround
        VINA_FOR(i, num_threads)
            create_thread(boost::ref(a));   // 创建工作线程，调用a.par->loop()
    }
    /**
     * @brief 同步并行计算，等待所有任务完成即可
     */
	void run(sz size_) {
		boost::mutex::scoped_lock self_lk(self);
        size = size_;
        finished = 0;
        started = 0;
        cond.notify_all(); // many things modified
        while(finished < size) // wait until processing of all elements is finished
            busy.wait(self_lk);
    }
    virtual ~parallel_for() {
        {
			boost::mutex::scoped_lock self_lk(self);
            destructing = true;
            cond.notify_all(); // destructing modified
        }
        join_all(); 
    }
private:
    /**
     * @brief 工作线程主循环（同步模式）,每个线程依次获取下一个可用的索引进行处理
     */
    void loop() {
        while(boost::optional<sz> i = get_next()) {
			(*m_f)(i.get());
			{
				boost::mutex::scoped_lock self_lk(self);
				++finished;         // 增加已完成任务计数
				busy.notify_one();  // 通知主线程
			}
		}
    }
    struct aux {
        parallel_for* par;
        aux() : par(NULL) {}
		void operator()() const { par->loop(); }
    };

    aux a;                          ///< 辅助对象实例
    const F* m_f;                   ///< 函数对象指针（parallel_mc_aux::operator()(parallel_mc_task& t)）
    boost::condition cond;          ///< 工作线程同步条件变量
    boost::condition busy;          ///< 主线程通知条件变量
    bool destructing;               ///< 析构标志
    sz size;                        ///< 任务总数
    sz started;                     ///< 已开始的任务数量
    sz finished;                    ///< 已完成的任务数量
    boost::mutex self;              ///< 互斥锁
    /**
     * @brief 获取下一个任务索引（线程安全）
     * @return 下一个任务索引的可选值，如果正在析构则返回空
     */
	boost::optional<sz> get_next() {
		boost::mutex::scoped_lock self_lk(self);
        while(!destructing && started >= size)
            cond.wait(self_lk);
		if(destructing) return boost::optional<sz>(); // NOTHING
		sz tmp = started;
        ++started;
        return tmp;
    }
};


/**
 * @brief 并行容器迭代器模板类
 * @tparam F 函数对象类型，用于处理容器元素
 * @tparam Container 容器类型
 * @tparam Input 输入类型（未使用）
 * @tparam Sync 同步模式标志，代码中只用到了true的情况，即同步模式
 */
template<typename F, typename Container, typename Input, bool Sync = false>
struct parallel_iter { 
	parallel_iter(const F* f, sz num_threads) : a(f), pf(&a, num_threads) {}
    /**
     * @brief 对容器进行并行迭代处理
     * @param v 待处理的容器引用
     */
	void run(Container& v) {
		a.v = &v;
		pf.run(v.size());
	}
private:
    /**
     * @brief 辅助结构体，将索引转换为容器元素访问
     */
	struct aux {
		const F* f;
		Container* v;
		aux(const F* f) : f(f), v(NULL) {}
        /**
         * @brief 处理指定索引的容器元素
         * @param i 元素索引
         */
		void operator()(sz i) const { 
			VINA_CHECK(v);
			(*f)((*v)[i]);  //相当于parallel_mc_aux中的void operator()(parallel_mc_task& t)，执行monte-carlo搜索
		}
	};
    aux a;                           ///< 辅助对象实例
    parallel_for<aux, Sync> pf;      ///< 并行for循环对象
};

#endif
