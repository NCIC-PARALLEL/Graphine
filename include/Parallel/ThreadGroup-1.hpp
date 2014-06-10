#ifndef _GRE_Sched_ThreadGroup_HPP
#define _GRE_Sched_ThreadGroup_HPP
/* -Master
 * -Thread Groups
 * --worker threads
 */
extern "C"{
	#define _GNU_SOURCE
	#include <sched.h>
}
#include <pthread.h>

#include "Config/SysConfig.hpp"
#include "Util/Locks.hpp"
#include "Util/vLocks.hpp"
#include "Util/FastForwardQueue.hpp"

#include <cstdlib>
#include <vector>
#include <string>

namespace GRE{
	//
	//Pre-declarision
	class ThreadGroup;
	class ThreadGroupShop;
	//
	//Thread Struct
	struct Thread {
		pthread_t tid;
		int cpuid;
		ThreadGroup& groupRef;

		Thread(ThreadGroup_t& _tg, int _cpuid)
			:tid(),groupRef(_tg), cpuid(_cpuid){}
	private:
		//forbid implicit construct.
		Thread(){}
	};
	typedef Thread Thread_t;
	//
	//
	class ThreadGroup {
	public:
		friend struct Thread;
		friend class ThreadGroupShop;
	private:
		typedef typename GRE::Util::vLock_AD vLock_t;
		typedef typename GRE::Util::FFQueue<byte*> FFQueue_t;
	public:
		ThreadGroup(const int gid, ThreadGroupShop& _global, const int _max_nThreads)
			:isInitialized(false), gid_(gid), global(_global), max_nThreads(_max_nThreads),nThreads(0),lock_() {
		}
		~ThreadGroup(){}
	private:
		//flag
		bool isInitialized;
		//Group ID
		int gid;
		//number of threads
		int nThreads;
		int max_nThreads;
		//Threads
		cpu_set_t cpuset;
		std::vector<Thread> threads;
		//sync structures
		vLock_t lock_;
		pthread_barrier_t barrier_;//32 bytes
		//Fast Channels among sockets
		std::vector<FFQueue*> inSocketChannels;
		std::vector<FFQueue*> outSocketChannels;
		//ref to global information structre
		ThreadGroupShop& global;
	public:
		ThreadGroupShop& getGlobal{
			return global;
		}
		inline int getGroupID() const
		{
			return gid;
		}
		inline int getNumThreads() const
		{
			return nThreads;
		}
		inline void barrier(){
			pthread_barrier_wait(&barrier_);
		}
		inline byte* getFromChannel(const int socketID){
			return inSocketChannels[socketID].dequeue();
		}
		inline bool putToChannel(const int socketID, const byte* buf){
			return outSocketChannels[socketID].enqueue(buf);
		}
		void setInSocketChannel(int src, FFQueue_t* channel){
			inSocketChannels[src] = channel;
		}
		void setOutSocketChannel(int dest, FFQueue_t* channel){
			outSocketChannels[dest] = channel;
		}
	public:
		void initialize(){
			const int nSockets = global.getNumSockets();

			//If use socket channels for communication, please set them later!
			inSocketChannels.resize(nSockets, NULL);
			outSocketChannels.resize(nSockets, NULL);

			CPU_SET(0, &cpuset);
			isInitialized = true;
		}
		void setAffinity(const cpu_set_t* _cpuset){
			//cpuset<--_cpuset
			memcpy(&cpuset, _cpuset, sizeof(cpu_set_t));
		}
		void launchThreads(pthread_attr_t* attr, void *(*fworker) (void *), const int n=0){
			assert(isInitialized == true);
			if(n!=0)
				nThreads = n;
			else
				nThreads = max_nThreads;

			pthread_t master = pthread_self();
			if(pthread_setaffinity_np(master, sizeof(cpu_set_t), &cpuset))
				perror("Fail to set affinity.");

			pthread_barrier_init(barrier_, NULL, nThreads);

			threads.reserve(n);
			for(int i=0; i<n; i++){
				Thread example(*this, i);
				threads.push_back(example);
				const int ret = pthread_create(&threads[i].tid, attr, fworker, &threads[i]);
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
			pthread_barrier_destroy(barrier_);
			threads.clear();
		}
		//After call finalize, the threadgroup can no longer be allowed to use. 
		void finalize(){
			CPU_ZERO(&cpuset);
			inSocketChannels.clear();
			outSocketChannels.clear();
			isInitialized = false;
		}
	};
	//
	//global control
	class ThreadGroupShop {
	public:
		typedef typename GRE::Util::FFQueue<byte*> FFQueue_t;
	private:
		//hardware information
		int nThreads;
		int nSockets;
		int nCoresPerSocket;
		int nThreadsPerCore;
		//infrastructure
		std::vector<ThreadGroup_t*> threadGroups;
		std::vector<FFQueue_t*> socketChannels;
		//sync instrument
		pthread_attr_t attr_;//56 bytes
		pthread_barrier_t barrier_;//32 bytes
		//
		bool isRunning;
		bool isInitialized;
	public:
		ThreadGroupShop():isRunning(false),isInitialized(false){}	
		~ThreadGroupShop(){}
	public:
		//Don't use the mannual config routine, if you don't exactly know the hardware parameters.
		void initialize(const int _nSocket = Sys_nSockets, const int _nCoresPerSocket = Sys_nCoresPerSocket,
			const int _nThreadsPerCore = Sys_nThreadsPerCore)
		{
			nSockets = _nSockets;
			nCoresPerSocket = _nCoresPerSocket;
			nThreadsPerCore = _nThreadsPerCore;
			nThreads = nSockets * nCoresPerSocket * nThreadsPerCore;
			//
			pthread_attr_init(&attr_);
			pthread_attr_setdetachstate(&attr_, PTHREAD_CREATE_JOINABLE);
			pthread_barrier_init(&barrier_,NULL,nThreads+1);	
			//
			socketChannels.resize(nSockets*nSockets, NULL);
			threadGroups.resize(nSockets, NULL);
			//create socket channels
			for(int i=0; i<nSockets; i++){
				////create infrastructure for each group
				for(int j=0; j<nSockets; j++){
					if(i != j)//no self to self
						socketChannels[i] = new FFQueue_t(NULL);
						socketChannels[i]->init();
				}
			}
			//create thread group infrastructure
			for(int i=0; i<nSockets; i++){
				const int nThreadsPerSocket = nThreadsPerCore*nCoresPerSocket;
				////Create and initialize threadgroup structure
				threadGroups[i] = new ThreadGroup(i, , *this, nThreadsPerSocket);	
				threadGroups[i].initialize();
				////Set socket channels
				for(int j=0; j<nSockets; j++){
					threadGroup[i].setInsocketChannel(j, socketChannels[i*nSockets+j]);//j-->i, parameter: src, channel
					threadGroup[i].setOutsocketChannel(j, socketChannels[j*nSokcets+i]);//i-->j, parameter: dest, channel
				}
			}
			isInitialized = true;
		}
		void launchThreads(void *(*fworker) (void *)){
			assert(isInitialized == true);
			cpu_set_t cpuset;
			const int step = nSockets*nCoresPerSocket;
			int idx = 0;
			for(int i=0; i<nSockets; i++){
				CPU_ZERO(&cpuset);
				for(int j=0; j<nThreadsPerCore; j++){
					idx = j*step + i*nCoresPerSocket;
					for(int j=0; j<nCoresPerSocket; j++){
						CPU_SET(idx,&cpuset);
						idx++;
					}
				}
				threadGroups[i].setAffinity(&cpuset);
				threadGroups[i].lauchThreads(&attr_, fworker);
				}
			}
			isRunning = true;
		}
		void joinThreads(){
			assert(isInitialized == true);
			//wait to finish
			if(isRunning){
				for(int i=0; i<nSockets; i++){
					threadGroup.joinThreads();
				}
				isRunning = false;
			}
		}
		void finalize(){
			assert(isInitialized == true);
			//if not join ,join.
			joinThreads();
			//destroy things
			for(int i=0; i<nSockets; i++){
				threadGroup.destroy();
			}
			//
			for(int i=0; i<nSockets; i++){
				////create infrastructure for each group
				for(int j=0; j<nSockets; j++){
					if(i != j)//no self to self
						socketChannels[i]->destroy();
						delete socketChannels[i];
				}
			}
			socketChannels.clear();
			//
			isInitialized = false;
		}
	};
	typedef ThreadGroupShop ProcessExt_t;
}
#endif
