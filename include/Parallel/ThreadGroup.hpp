#ifndef _GRE_ThreadGroup_HPP
#define _GRE_ThreadGroup_HPP
/* -Master
 * -Thread Groups
 * --worker threads
 */
extern "C"{
	#ifndef _GNU_SOURCE
	#define _GNU_SOURCE
	#endif
	#include <sched.h>
	#include <string.h>
}
#include <pthread.h>

#include "Config/SysConfig.hpp"
#include "Util/Locks.hpp"
#include "Util/vLock.hpp"
#include "Util/FastForwardQueue.hpp"

#include <cstdlib>
#include <vector>
#include <cassert>

namespace GRE{
	//
	//Pre-declarision
	//
	class ThreadGroupShop;
	//
	class ThreadGroup{
	private:
		struct Thread{
			pthread_t tid;
			int id;
			ThreadGroup* ptrThreadGroup;
			void* rtm;
			//
			Thread(ThreadGroup* _tg, int _id, void* _rtm=NULL)
				:tid(),ptrThreadGroup(_tg), id(_id), rtm(_rtm){}
		private:
			//forbid implicit construct.
			Thread(){}
		public:
			//obsolete
			inline ThreadGroup& getRefThreadGroup(){
				return (*ptrThreadGroup);
			}
			//obsolete
			inline ThreadGroupShop& getRefThreadGroupShop(){
				return *(ptrThreadGroup->global);
			}
			//recommented interface replacing getRefxx
			inline ThreadGroup& refThreadGroup(){
				return (*ptrThreadGroup);
			}
			//recommented interface replacing getRefxx
			inline ThreadGroupShop& refThreadGroupShop(){
				return *(ptrThreadGroup->global);
			}
		};
	public:
		friend class ThreadGroupShop;
	public:
		typedef Thread Thread_t;
		typedef GRE::Util::vLock_ID vLock_t;
		typedef GRE::Util::FFQueue<byte*> FFQueue_t;
	public:
		ThreadGroup(ThreadGroupShop* _global, const int _gid, const int _max_nThreads)
			:isInitialized(false),global(_global),gid(_gid),max_nThreads(_max_nThreads),nThreads(0),
			vlock_(),inSocketChannel(NULL),outSocketChannels(){}
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
		vLock_t vlock_;
		pthread_barrier_t barrier_;//32 bytes
		//Fast Channels among sockets
		FFQueue_t* inSocketChannel;
		std::vector<FFQueue_t*> outSocketChannels;
		//ref to global information structre
		ThreadGroupShop* global;
	public:
		//obsolete
		inline ThreadGroupShop& getRefGlobal(){
			return (*global);
		}
		//use this one
		inline ThreadGroupShop& refGlobal(){
			return (*global);
		}
		//obsolete
		inline vLock_t& getRefvLock(){
			return vlock_;
		}
		//use this one
		inline vLock_t& refvLock(){
			return vlock_;
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
		inline bool isLeader(int tid){
			return (tid%nThreads==0);
		}
		//
		inline byte* getFromChannel(){
			return inSocketChannel->dequeue();
		}
		inline bool putToChannel(const int socketID, byte* buf){
			return outSocketChannels[socketID]->enqueue(buf);
		}

		void setInSocketChannel(FFQueue_t* channel){
			inSocketChannel = channel;
		}
		void setOutSocketChannel(const int dest, FFQueue_t* channel){
			if(dest < outSocketChannels.size())
				outSocketChannels[dest] = channel;
			else//maybe not correct! it depends.
				outSocketChannels.push_back(channel);
		}
	public:
		void initialize(){
			CPU_SET(0, &cpuset);
			inSocketChannel = NULL;
			isInitialized = true;
		}
		void setAffinity(const cpu_set_t* _cpuset){
			//cpuset<--_cpuset
			memcpy(&cpuset, _cpuset, sizeof(cpu_set_t));
		}
		void launchThreads(pthread_attr_t* attr, void *(*fworker) (void *), void* rtm, const int n=0){
			if(!isInitialized)
				initialize();
			assert(isInitialized == true);
			if(n>0)
				nThreads = n;
			else
				nThreads = max_nThreads;

			pthread_t master = pthread_self();
			if(pthread_setaffinity_np(master, sizeof(cpu_set_t), &cpuset))
				std::cerr<<"Fail to set affinity."<<std::endl;

			pthread_barrier_init(&barrier_, NULL, nThreads);

			threads.reserve(nThreads);
			for(int i=0; i<nThreads; i++){
				threads.push_back(Thread(this, nThreads*gid+i, rtm));
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
			pthread_barrier_destroy(&barrier_);
			threads.clear();
		}
		//After call finalize, the threadgroup can no longer be allowed to use. 
		void finalize(){
			CPU_ZERO(&cpuset);
			inSocketChannel = NULL;
			outSocketChannels.clear();
			isInitialized = false;
		}
	};
	//
	//global control
	class ThreadGroupShop {
	public:
		typedef ThreadGroup ThreadGroup_t;
		typedef GRE::Util::FFQueue<byte*> FFQueue_t;
		typedef GRE::Util::vLock_ID vLock_t;
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
		//
		vLock_t vlock_;
	public:
		ThreadGroupShop():isRunning(false),isInitialized(false),vlock_(){}	
		~ThreadGroupShop(){}
	public:
		//Don't use the mannual config routine, if you don't exactly know the hardware parameters.
		void initialize(const int _nSockets = Sys_nSockets, const int _nCoresPerSocket = Sys_nCoresPerSocket,
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
			socketChannels.resize(nSockets, NULL);
			threadGroups.resize(nSockets, NULL);
			//create socket channels
			for(int i=0; i<nSockets; i++){
				socketChannels[i] = new FFQueue_t(NULL);
				socketChannels[i]->init(2048);
			}
			//create thread group infrastructure
			const int nThreadsPerSocket = nThreadsPerCore*nCoresPerSocket;
			for(int i=0; i<nSockets; i++){
				////Create and initialize threadgroup structure
				threadGroups[i] = new ThreadGroup(this, i, nThreadsPerSocket);	
				threadGroups[i]->initialize();
				////Set socket channels
				threadGroups[i]->setInSocketChannel(socketChannels[i]);//j-->i, parameter: src, channel
				for(int j=0; j<nSockets; j++){
					threadGroups[i]->setOutSocketChannel(j, socketChannels[j]);//i-->j, parameter: dest, channel
				}
			}
			isInitialized = true;
		}
		void launchThreads(void *(*fworker) (void *), void* rtm=NULL){
			if(!isInitialized)
				initialize();//with default parameters
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
				threadGroups[i]->setAffinity(&cpuset);
				threadGroups[i]->launchThreads(&attr_, fworker, rtm);
			}
			isRunning = true;
		}
		void joinThreads(){
			assert(isInitialized == true);
			//wait to finish
			if(isRunning){
				for(int i=0; i<nSockets; i++){
					threadGroups[i]->joinThreads();
				}
				isRunning = false;
			}
		}
		void finalize(){
			if(isInitialized){
				//if not join ,join.
				joinThreads();
				//destroy things
				for(int i=0; i<nSockets; i++){
					threadGroups[i]->finalize();
					delete threadGroups[i];
				}
				threadGroups.clear();
				//
				for(int i=0; i<nSockets; i++){
					socketChannels[i]->destroy();
					delete socketChannels[i];
				}
				socketChannels.clear();
				//
				pthread_barrier_destroy(&barrier_);	
				//
				isInitialized = false;
			}
		}
		//Helper
		inline int getNumSockets(){
			return nSockets;
		}
		inline int getNumGroups(){
			return nSockets;
		}
		inline int getNumThreads(){
			return nThreads;
		}
		inline void barrier(){
			pthread_barrier_wait(&barrier_);
		}
		//obsolete
		inline vLock_t& getRefvLock(){
			return vlock_;
		}
		//use this one
		inline vLock_t& refvLock(){
			return vlock_;
		}
	};
	typedef ThreadGroupShop ProcessExt_t;
}
#endif
