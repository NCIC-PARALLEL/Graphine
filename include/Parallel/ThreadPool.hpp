#ifndef _GRE_ThreadPool_HPP
#define _GRE_ThreadPool_HPP
/* -Master
 * -Thread Pools
 * --worker threads
 */
#include <pthread.h>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "Config/SysConfig.hpp"
#include "Util/Locks.hpp"
#include "Util/vLock.hpp"

#include <cassert>
namespace GRE{
	//
	//
	//
	class ThreadPool{
	private:
		//Thread Struct
		struct Thread{
			pthread_t tid;
			int id;
			ThreadPool* ptrThreadPool;
			void* rtm;
			//
			Thread(ThreadPool* _tp, int _id, void* _rtm=NULL)
				:tid(), ptrThreadPool(_tp), id(_id), rtm(_rtm){}
			//
			//Helper ROutines
			inline ThreadPool& refThreadPool(){
				return (*ptrThreadPool);
			}
		};
	public:
		typedef Thread Thread_t;
		typedef GRE::Util::vLock_ID vLock_t;
	public:
		ThreadPool(const int _max_nThreads)
			:isInitialized(false), max_nThreads(_max_nThreads), nThreads(0), lock_() {
		}
		~ThreadPool(){}
	private:
		//flag
		bool isInitialized;
		//number of threads
		int nThreads;
		int max_nThreads;
		//Threads
		std::vector<Thread_t> threads;
		//sync structures
		vLock_t lock_;
		pthread_attr_t attr_;//56 bytes
		pthread_barrier_t barrier_;//32 bytes
	public:
		inline int getNumThreads() const
		{
			return nThreads;
		}
		inline void barrier(){
			pthread_barrier_wait(&barrier_);
		}
		inline vLock_t& getRefvLock(){
			return lock_;
		}
		inline vLock_t& refvLock(){
			return lock_;
		}

	public:
		void initialize(){
			pthread_attr_init(&attr_);
			pthread_attr_setdetachstate(&attr_, PTHREAD_CREATE_JOINABLE);
			isInitialized = true;
		}
		void launchThreads(void *(*fworker) (void *), void* data, const int n=0){
			if(!isInitialized)
				initialize();
			if(n > 0)
				nThreads = n;
			else
				nThreads = max_nThreads;

			pthread_barrier_init(&barrier_, NULL, nThreads+1);

			threads.reserve(nThreads);
			for(int i=0; i<nThreads; i++){
				threads.push_back(Thread(this, i, data));
				const int ret = pthread_create(&threads[i].tid, &attr_, fworker, &threads[i]);
				if(ret){
					std::cerr<<"Err: fail to create a thread."<<std::endl;
				}
			}
		}
		void joinThreads(){
			assert(isInitialized == true);
			assert(nThreads > 0);

			for(int i=0; i<nThreads; i++){
				const int ret = pthread_join(threads[i].tid, NULL);
				if(ret){
					std::cerr<<"Err: fail to join a thread."<<std::endl;
				}
			}
			nThreads = 0;
			pthread_barrier_destroy(&barrier_);
			threads.clear();
		}
		//After call finalize, the threadgroup can no longer be allowed to use. 
		void finalize(){
			pthread_attr_destroy(&attr_);
			isInitialized = false;
		}
	};
}
#endif
