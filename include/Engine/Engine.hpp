#ifndef _GRE_Engine_HPP
#define _GRE_Engine_HPP

#include "Types.hpp"
#include "Status.hpp"
#include "Message.hpp"
#include "Config/UsrConfig.hpp"

#include "Parallel/ThreadPool.hpp"
#include "Parallel/ThreadGroup.hpp"
#include "Parallel/ThreadContext.hpp"

#include "Graph/DistributedGraph.hpp"
#include "Graph/VertexSet.hpp"

#include "Util/Util.hpp"
#include "Util/Buffer.hpp"
#include "Util/BufferPool.hpp"
#include "Util/vLock.hpp"
#include "Util/Timer.hpp"

#include "Comm/BufferingComm.hpp"
#include "Comm/CommUnderlying.hpp"

#include <vector>
#include <iostream>
#include <sched.h>

namespace util = GRE::Util;

#define vbufSize 128
#define scheduleGrain 8

namespace GRE{
	typedef enum ParallelMode{
		Par_ThreadPool,
		Par_ThreadGroup,
		Par_Other
	} ParallelMode_t;
	typedef enum CtrlState{
		PS_Waiting,
		PS_Running,
		PS_Idle,
		PS_Done,
		PS_InitAll,
		PS_InitSet,
		PS_Phase_Sync_Scatters,
		PS_Phase_Sync_Combiners,
		PS_Phase_Compute,
		PS_Phase_Apply
	} CtrlState_t;
	//BSP-mode
	template <typename VertexCode>
	class SynchronousEngine{
	private:
		typedef typename VertexCode::ScatterData_t ScatterData_t;
		typedef typename VertexCode::CombineData_t CombineData_t;
		typedef GRE::Msg<ScatterData_t, vid_t> ScatterMsg_t;
		typedef GRE::Msg<CombineData_t, vid_t> CombineMsg_t;
		typedef GRE::Msg<CombineData_t, lvid_t> SocketMsg_t;

		typedef GRE_Status Status_t;
		typedef GRE::ThreadGroupShop ThreadGroupShop_t;
		typedef GRE::ThreadGroup ThreadGroup_t;
		typedef GRE::ThreadPool ThreadPool_t;

		typedef GRE::Graph::VertexSet VertexSet_t;
		typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;
		typedef GRE::Graph::GraphInCSR LocalGraph_t;
		typedef GRE::Graph::DistributedGraph<LocalGraph_t> DistributedGraph_t;

		typedef GRE::Util::Bitmap Bitmap_t;
		typedef GRE::Util::FFQueue<byte*> FFQueue_t;
		typedef GRE::Util::BufferPool BufferPool_t;
		typedef GRE::Util::parBufferPool parBufferPool_t;
		typedef GRE::Util::FormatedBuffer<GRE::COMM::Protocol::Header, ScatterMsg_t> ScatterFormatedBuffer_t;
		typedef GRE::Util::FormatedBuffer<GRE::COMM::Protocol::Header, CombineMsg_t> CombineFormatedBuffer_t;


		typedef GRE::COMM::Underlying CommUnderlying_t;
		typedef GRE::COMM::BufComm<parBufferPool_t> BufComm_t;

		typedef GRE::SynchronousEngine<VertexCode> SELF_t;
		typedef GRE::ThreadContext<ThreadPool_t, VertexSet_t, SELF_t> TP_Context_t;
		typedef GRE::TG_ThreadContext<VertexSet_t, SELF_t> TG_Context_t;

		typedef GRE::TSCTimer Timer_t;
	public:
		typedef GRE::Util::FormatedBuffer<GRE::COMM::Protocol::Header, SocketMsg_t> SocketFormatedBuffer_t;
	private:
		//Don't allow default construction without initialization.
		SynchronousEngine(){}
	public:
		SynchronousEngine(DistributedGraph_t& _dg, CommUnderlying_t& _comm, BufferPool_t& _pool, VertexCode& _vertexCode)
			:graph(_dg), vertexCode(_vertexCode), currSet(), nextSet(),incomingQueue(NULL), outgoingQueue(NULL), dc(_comm, interBufPool)
		{
			intraBufPool.inherit(_pool, 8192*2);
			interBufPool.inherit(_pool, 4096);

			incomingQueue.init(4096);
			outgoingQueue.init(4096);

			currSet.initialize(_dg, true);
			nextSet.initialize(_dg, false);

			isInitialized = false;

			max_nSupersteps = 1024;
			syncState = PS_Waiting;

			nFinishedGroups = 0;
		}

		void reset(){
			superstep = 0;
			nFinishedGroups = 0;
			syncState = PS_Waiting;
			incomingQueue.reset();
			outgoingQueue.reset();
			currSet.clear();
			nextSet.clear();
			interBufPool.reset();
			intraBufPool.reset();
			vertexCode.free();
			vertexCode.alloc();

			initVertexTimer.reset(); 
			syncScatterTimer.reset(); 
			syncCombinerTimer.reset(); 
			localComputeTimer.reset(); 
			localApplyTimer.reset(); 
		}
		~SynchronousEngine(){
			incomingQueue.destroy();
			outgoingQueue.destroy();

			currSet.destroy();
			nextSet.destroy();
	
			intraBufPool.destroy();
			interBufPool.destroy();

			isInitialized = false;
		}
		//
		//Configuration
		//
		void configMaxSupersteps(const long n=0){
			if(n>0)
				max_nSupersteps = n;
			else
				max_nSupersteps = 1024;
		}
		void initialize(){
			if(!isInitialized){
				isInitialized = true;
				vertexCode.alloc();
			}

			initVertexTimer.init(); 
			syncScatterTimer.init(); 
			syncCombinerTimer.init(); 
			localComputeTimer.init(); 
			localApplyTimer.init(); 
		}
		void finalize(){
			if(isInitialized){
				isInitialized = false;
				vertexCode.free();
			}
		}
		//
		//
		//Core Computation
		//We offer two kinds of parallel methods, i.e. thread pool and thread groups.
		//
		void run(SourceVertexSet_t& srcSet, const ParallelMode parMode=Par_ThreadPool){
			run<Out>(srcSet, parMode);
		}
		template<Dir dir>
		void run(SourceVertexSet_t& srcSet, const ParallelMode parMode=Par_ThreadPool){
			if(superstep!=0){
				superstep = 0;
				nFinishedGroups = 0;
				syncState = PS_Waiting;
				incomingQueue.reset();
				outgoingQueue.reset();
				currSet.clear();
				nextSet.clear();
				interBufPool.reset();
				intraBufPool.reset();
			}
			sourceVertexSet_ptr = &srcSet;
			if(parMode==Par_ThreadPool)
				run_tp<dir>();
			else if(parMode==Par_ThreadGroup)
				run_tg<dir>();
			else {
				//should not be here!
				std::cerr<<"Fetal Err: this kind of multi-threading mode is Not Implemented."<<std::endl;
			}
		}
		//
		//ThreadPool-based Parallel computation
		//
		template<Dir dir>
		Status_t run_tp(){
			ThreadPool_t lc(GRE::MAX_NUM_THREADS);
			lc.initialize();
			//run worker
			lc.launchThreads(run_worker_tp<dir>, static_cast<void*>(this), 12);
			//run master
			Status_t retMaster = run_master(lc);
			//parse the status of master thread
			if(retMaster != GRE_SUCCESS){
				GRE_Status_Parser parser;
				std::cerr<<"Master return with "<<parser.parse(retMaster)<<std::endl;
			}
			//join all threads in pool
			lc.joinThreads();
			//finalize
			lc.finalize();
		}
		//
		//ThreadGroup-based Parallel computation
		//
		template<Dir dir>
		Status_t run_tg(){
			ThreadGroupShop_t lc;
			lc.initialize(2, 6, 1);
			//run workers
			lc.launchThreads(run_worker_tg<dir>, static_cast<void*>(this));
			//run master
			Status_t retMaster = run_master(lc);
			//parse the status of master thread
			if(retMaster != GRE_SUCCESS){
				GRE_Status_Parser parser;
				std::cerr<<"Master return with "<<parser.parse(retMaster)<<std::endl;
			}
			//join all threads in pool
			lc.joinThreads();
			//finalize
			lc.finalize();
		}
		//template<Dir dir>
		template <typename LocalControl_t>
		Status_t run_master(LocalControl_t& lc) {//I am master.
			const int nProcs = dc.nProcs();
			const int procID = dc.procID();
			////init
			SourceVertexSet_t& srcSet = *(sourceVertexSet_ptr);
			if(srcSet.isALL()){//full set
				syncState = PS_InitAll;
			} else {//subset, should have small amount of vertices
				currSet.setAndInsert(srcSet.vertices);
				syncState = PS_InitSet;
			}
			lc.barrier();//-------------------------------------------2
			initVertexTimer.start();
			//wait workers to finish...............
			lc.barrier();//-----------------------------------------------3
			initVertexTimer.stop();
			
			////level-iterations
			currSet.refresh();
			dc.barrier();//----------------------------------------------g1
			//std::cerr<<procID<<": Number of ready vertices..."<<currSet.size()<<std::endl;
			//std::cout<<procID<<":"<<graph.get_num_masters()<<"::"<<graph.get_num_scatters()<<"::"<<graph.get_num_combiners()<<std::endl;
			for(superstep = 0; superstep < max_nSupersteps; superstep++){
				//std::cerr<<superstep<<": Number of ready vertices..."<<currSet.size()<<std::endl;
				syncState = PS_Phase_Sync_Scatters;
				lc.barrier();//----------------------------------------------4-0
				///////////////////////////Phase-1: sync scatters////////////////////////////////
				syncScatterTimer.start();
				do {
					int nSendActive = nProcs - 1;
					int nRecvActive = nProcs - 1;
					dc.launchAllRecv();
					while(nRecvActive > 0 || nSendActive > 0) {
						//send
						while(nSendActive > 0){
							byte* addr = outgoingQueue.dequeue_serial();
							if(addr){
								ScatterFormatedBuffer_t fbuf(addr, interBufPool.bufSize());
								if(fbuf.getHeader()->isEnd())
									nSendActive--;
								const int dest = fbuf.getHeader()->getExtra();
								fbuf.getHeader()->setExtra(procID);
								dc.isend(fbuf, dest);
							} else 
								break;
						}
						//recv
						if(dc.checkAllRecv()){
							dc.checkAllSend();
							typename BufComm_t::incomingQueue_t& queue = dc.getIncomingQueue();
							byte* addr;
							while(queue.get(addr)){
								ScatterFormatedBuffer_t fbuf(addr, interBufPool.bufSize());
								if(fbuf.getHeader()->isEnd()){
									dc.cancelRecv(fbuf.getHeader()->getExtra());
									nRecvActive--;
								}
								while(!incomingQueue.enqueue_serial(addr));//should never fails!
							}
						} 
					}	
					while(dc.checkAllSend());
					incomingQueue.waitForEmpty();
				} while(0);

				syncState = PS_Phase_Compute;
				dc.barrier();//----------------------------------------------g2
				lc.barrier();//-------------------------------------------4-2
				syncScatterTimer.stop();

				//>>nextSet
				///////////////////////////Phase-2: Main Computation////////////////////////////////
				//std::cerr<<superstep<<":"<<procID<<" | Number of scatters..."<<currSet.num_new()<<std::endl;
				//std::cerr<<superstep<<":"<<procID<<" | Number of vertices..."<<currSet.num_ready()<<std::endl;
				currSet.refresh();
				nextSet.clear();
				lc.barrier();//-------------------------------------------4-3
				localComputeTimer.start();
				//compute on local vertices
				syncState = PS_Phase_Sync_Combiners;
				lc.barrier();//-------------------------------------------4-4
				localComputeTimer.stop();
				//>>nextSet
				///////////////////////////Phase-3: Sync Combiners////////////////////////////////
				//

				syncCombinerTimer.start();
				do {
					int nSendActive = nProcs - 1;
					int nRecvActive = nProcs - 1;
					dc.launchAllRecv();
					while(nRecvActive > 0 || nSendActive >0){
						//send
						if(nSendActive > 0){
							byte* addr;
							while(outgoingQueue.dequeue_serial(addr)){
								CombineFormatedBuffer_t fbuf(addr, interBufPool.bufSize());
								if(fbuf.getHeader()->isEnd())
									nSendActive--;
								const int dest = fbuf.getHeader()->getExtra();
								fbuf.getHeader()->setExtra(procID);
								dc.isend(fbuf, dest);
							}
						}
						//recv
						if(nRecvActive > 0){
							dc.checkAllSend();
							if(dc.checkAllRecv()){
								typename BufComm_t::incomingQueue_t& queue = dc.getIncomingQueue();
								byte* addr;
								while(queue.get(addr)){
									CombineFormatedBuffer_t fbuf(addr, interBufPool.bufSize());
									if(fbuf.getHeader()->isEnd()){
										dc.cancelRecv(fbuf.getHeader()->getExtra());
										nRecvActive--;
									}
									while(!incomingQueue.enqueue_serial(addr))
										;//should never fails!
								}
							} 
						}
					}	
					while(dc.checkAllSend());
					incomingQueue.waitForEmpty();
				} while(0);
				syncState = PS_Phase_Apply;
				lc.barrier();//------------------------------------------4-5
				syncCombinerTimer.stop();
				//>>currSet
				///////////////////////////Phase-4: Apply////////////////////////////////
				currSet.clear_set();

				lc.barrier();//------------------------------------------4-6
				localApplyTimer.start(); 
				//compute on local vertices
				lc.barrier();//------------------------------------------4-7
				localApplyTimer.stop(); 

				//5:vote to halt
				currSet.refresh();
				//std::cerr<<superstep<<"|"<<procID<<":Number of ready vertices..."<<currSet.size()<<std::endl;
				//not necessary to sync
				//dc.barrier();//----------------------------------------------extra
				//
				vertexCode.reduce(dc);
				long myStat = currSet.size();
				long allStat = 0;
				//implicit global synshronization
				dc.refUnderlying().allReduce(&myStat, &allStat, 1, MPI_LONG, MPI_SUM);
				if(procID == 0){
					std::cout<<"Superstep-"<<superstep<<": Number of active vertices..."<<allStat<<std::endl;
				}
				if(allStat == 0)
					break;
			}//end for
			syncState = PS_Done;
			lc.barrier();//----------------------------------------------5(4-1)
			output_statistics(dc);
			return GRE_SUCCESS;
		}
		//
		//TP-Workers
		//Driven by master thread.
		//
		template<Dir dir>
		static void* run_worker_tp(void* arg){
			ThreadPool::Thread_t* thread = static_cast<ThreadPool::Thread_t*>(arg);
			Status_t ret = (static_cast<SELF_t*>(thread->rtm))->fn_worker_tp<dir>(*thread);
			if(ret!=GRE_SUCCESS){
				GRE_Status_Parser parser;
				std::cerr << "Err: workers(thread"<<thread->id
					<<") failure (failure number: "<<parser.parse(ret)<<")."<<std::endl;
			}
			pthread_exit(0);
		}
		template<Dir dir>
		Status_t fn_worker_tp(ThreadPool_t::Thread_t& thread){//We are workers!
			ThreadPool_t& lc = thread.refThreadPool();
			const int tid = thread.id;
			const int nThreads = lc.getNumThreads();
			const int nProcs = dc.nProcs();
			const int procID = dc.procID();
			LocalGraph_t& localGraph = graph.refLocalGraph();
			
			//construct thread context
			TP_Context_t ctx(lc, nextSet, this);
			
			lc.barrier();//------------------------------------2
			//initilization
			//do in parallel, statically divide the work(vertices).
			if(syncState == PS_InitAll){
				const int chunkSize = graph.get_num_masters()/nThreads + (graph.get_num_masters()%nThreads==0? 0:1);
				const lvid_t begin = tid*chunkSize;
				const lvid_t end = (tid+1)*chunkSize < graph.get_num_masters()? (tid+1)*chunkSize : graph.get_num_masters();
				lvid_t vbuf[vbufSize];
				int nv=0;
				for(lvid_t lv=begin; lv<end; lv++){
					if(vertexCode.init(lv)){
						currSet.set(lv);
						vbuf[nv++]=lv;
						if(nv==vbufSize){
							currSet.insert(vbuf, nv);
							nv = 0;
						}
					}
				}
				if(nv > 0){
					currSet.insert(vbuf, nv);
					nv = 0;
				}
			}
			//or
			if(syncState == PS_InitSet && tid==0){
				SourceVertexSet_t& srcSet = *(sourceVertexSet_ptr);
				const size_t n = srcSet.vertices.size();
				for(int i=0; i<n; i++){
					const vid_t v = srcSet.vertices[i];
					if(graph.get_owner(v)==procID){
						const lvid_t lv = graph.get_lvid_master(v);
						if(vertexCode.init(lv)){
							currSet.setAndInsert(lv);
						}
					}
				}
			}
			lc.barrier();//----------------------------------------3

			//Main-loop
			while(1){
				lc.barrier();//----------------------------------------------4-0
				if(syncState == PS_Done) break;
				//1:Sync scatters
				//update native scatters
				///////////////////////////Phase-1////////////////////////////////
				//assert(syncState == PS_Phase_Sync_Scatters);
				//1: send
				const int nThreads1=nThreads-1;
				if(tid!=nThreads1)
				for(int p = tid; p < nProcs; p+=nThreads1){
					if(p==procID) continue;
					ScatterFormatedBuffer_t fbuf;
					byte* addr;
					do {
						addr = recvScatters(incomingQueue,currSet);
						if(!addr)
							addr = interBufPool.get_address();
					} while(!addr);//until get one. should not fail.
					fbuf.format(addr, interBufPool.bufSize());

					std::vector<lvid_t>& activeSet = currSet.exportSet(); 
					const int64_t beginIdx = 0;
					const int64_t endIdx = currSet.num_ready();
					if(dir==Out){
						Bitmap_t& scatterMap = graph.scatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(scatterMap.test(lv)){//lv has scatters in remote machine p				
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					} 
					if(dir==In){
						Bitmap_t& rscatterMap = graph.rscatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(rscatterMap.test(lv)){//lv has rscatters(combiners) in remote machine p				
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==InOut){
						Bitmap_t& scatterMap = graph.scatterMaps[p];
						Bitmap_t& rscatterMap = graph.rscatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(scatterMap.test(lv)||rscatterMap.test(lv)){
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					fbuf.getHeader()->setExtra(p);
					fbuf.getHeader()->setEnd();
					while(!outgoingQueue.enqueue(fbuf.address()))
						;//never fails! to ensure this.
				}
				//2:recv
				if(tid!=nThreads1)
				while(syncState == PS_Phase_Sync_Scatters){
					byte* addr = recvScatters(incomingQueue,currSet);
					if(addr){
						interBufPool.put_address(addr);
					} else 
						sched_yield();		
				}
				lc.barrier();//------------------------------------4-2
				//wait master to update currSet.
				//2:Compute
				///////////////////////////Phase-2: main Computation////////////////////////////////
			//	assert(syncState == PS_Phase_Compute);
				lc.barrier();//------------------------------------4-3
				do {
					lvid_t buf[scheduleGrain];
					while(1){
						const int n = currSet.fetch(buf, scheduleGrain);
						if(n > 0){
							if(dir==Out || dir==InOut)
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								LocalGraph_t::edge_iterator it = localGraph.edge_begin(lu);
								const LocalGraph_t::edge_iterator end = localGraph.edge_end(lu);
								for(;it!=end; it++){
									const lvid_t lv = *it;
									vertexCode.scatter(ctx, lu, lv, it);
								}
							}
							if(dir==In || dir==InOut)
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								LocalGraph_t::redge_iterator it = localGraph.redge_begin(lu);
								const LocalGraph_t::redge_iterator end = localGraph.redge_end(lu);
								for(;it!=end; it++){
									const lvid_t lv = *it;
									vertexCode.scatter(ctx, lu, lv, it);
								}
							}
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								if(vertexCode.assert_to_halt(lu))
									currSet.unset(lu);
							}
						} else
							break;
					}
				} while(0);
				lc.barrier();//------------------------------------4-4
				//3:Sync combiners
				///////////////////////////Phase-3: sync combiners////////////////////////////////
			//	assert(syncState == PS_Phase_Sync_Combiners);
				if(tid!=nThreads1)
				for(int p = tid; p < nProcs; p+=nThreads1){
					if(p==procID) continue;
					CombineFormatedBuffer_t fbuf;
					byte* addr;
					do {
						addr = recvCombiners(ctx, incomingQueue, nextSet);
						if(!addr)
							addr = interBufPool.get_address();
					} while(!addr);//until get one. should not fail.
					fbuf.format(addr, interBufPool.bufSize());


					if(dir==Out){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& combinerMap = graph.combinerMaps[p];
						Bitmap_t::wordIterator it1 = combinerMap.begin();
						//const Bitmap_t::wordIterator end = activeBitset.getWordIterator(graph.get_end_lvid_combiner());
						//const Bitmap_t::wordIterator end1 = combinerMap.end();
						const lvid_t maxID = graph.get_end_lvid_combiner();
						for(lvid_t base = graph.get_begin_lvid_combiner();base < maxID; it++, it1++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==In){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& rcombinerMap = graph.rcombinerMaps[p];
						Bitmap_t::wordIterator it1 = rcombinerMap.begin();
						const lvid_t maxID = graph.get_end_lvid_scatter();
						for(lvid_t base = graph.get_begin_lvid_scatter();base < maxID; it++, it1++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==InOut){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& combinerMap = graph.combinerMaps[p];
						Bitmap_t::wordIterator it1 = combinerMap.begin();
						Bitmap_t& rcombinerMap = graph.rcombinerMaps[p];
						Bitmap_t::wordIterator it2 = rcombinerMap.begin();
						const lvid_t maxID = graph.get_end_lvid_combiner();
						for(lvid_t base = graph.get_begin_lvid_combiner();base < maxID; it++, it1++, it2++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1 | *it2), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					
					fbuf.getHeader()->setExtra(p);
					fbuf.getHeader()->setEnd();
					while(!outgoingQueue.enqueue(fbuf.address()))
						;//never fails! to ensure this.
				}

				//2:recv
				if(tid!=nThreads1)
				while(syncState == PS_Phase_Sync_Combiners){
					byte* addr = recvCombiners(ctx, incomingQueue,nextSet);
					if(addr){
						interBufPool.put_address(addr);
					} else 
						sched_yield();				
				}
				lc.barrier();//------------------------------------------4-5
				///////////////////////////Phase-4////////////////////////////////
				lc.barrier();//------------------------------------------4-6
				assert(syncState == PS_Phase_Apply);
				const size_t masterIdx_begin = graph.get_begin_lvid_master();
				const size_t masterIdx_end = graph.get_end_lvid_master();
				Bitmap_t& activeBitset = nextSet.exportBitset(); 
				const size_t nWords = util::align(masterIdx_end - masterIdx_begin)/64;
				//apply
				for(size_t i = masterIdx_begin/64+tid; i<nWords; i+=nThreads){
					uint64_t tmp = activeBitset.getWord(i);
					lvid_t lv = i*64;
					while(tmp){
						if(tmp & 0x1UL){
							if(vertexCode.apply(lv)) currSet.set(lv);
						}
						lv++;
						tmp>>=1;
					}
				}
				//refresh currSet for next superstep
				lvid_t vbuf[vbufSize];
				int nv = 0;
				Bitmap_t& activeBitset1 = currSet.exportBitset(); 
				for(size_t i = masterIdx_begin/64+tid; i<nWords; i+=nThreads){
					uint64_t tmp = activeBitset1.getWord(i);
					lvid_t lv = i*64;
					while(tmp){
						if(tmp & 0x1UL){
							vbuf[nv++] = lv;
							if(nv == vbufSize){
								currSet.insert(vbuf, nv);
								nv = 0;
							}
						}
						lv++;
						tmp>>=1;
					}
				}
				if(nv > 0){
					currSet.insert(vbuf, nv);
					nv = 0;
				}
				lc.barrier();//------------------------------------------4-7
				/*
				if(!vertexCode.aggregator){
					VertexCode::Aggregator agg;
					const int chunkSize = graph.get_num_masters()/nThreads + (graph.get_num_masters()%nThreads==0? 0:1);
					const lvid_t begin = tid*chunkSize;
					const lvid_t end = (tid+1)*chunkSize < graph.get_num_masters()? (tid+1)*chunkSize : graph.get_num_masters();
					for(lvid_t lv=begin; lv<end; lv++){
						agg.aggregate(lv);
					}
					lc.barrier();//------------------------------------------4-8
				}
				*/
			}
			return GRE_SUCCESS;
		}

		//
		//TG-Workers
		//Driven by master thread.
		//
		template<Dir dir>
		static void* run_worker_tg(void* arg){
			ThreadGroup_t::Thread_t* thread = static_cast<ThreadGroup_t::Thread_t*>(arg);
			Status_t ret = (static_cast<SELF_t*>(thread->rtm))->fn_worker_tg<dir>(*thread);
			if(ret!=GRE_SUCCESS){
				GRE_Status_Parser parser;
				std::cerr << "Err: workers(thread"<<thread->id
					<<") failure (failure number: "<<parser.parse(ret)<<")."<<std::endl;
			}
			pthread_exit(0);
		}
		template<Dir dir>
		Status_t fn_worker_tg(ThreadGroup_t::Thread_t& thread){//We are workers!
			ThreadGroupShop_t& lc = thread.refThreadGroupShop();
			ThreadGroup_t& group = thread.refThreadGroup();
			const int tid = thread.id;
			const int nThreads = lc.getNumThreads();
			const int nGroups = lc.getNumGroups();
			const int nProcs = dc.nProcs();
			const int procID = dc.procID();
			LocalGraph_t& localGraph = graph.refLocalGraph();
		
			//construct thread context
			TG_Context_t ctx(thread, nextSet, this);
			
			lc.barrier();//------------------------------------2
			//initilization
			//do in parallel, statically divide the work(vertices).
			if(syncState == PS_InitAll){
				const int chunkSize = graph.get_num_masters()/nThreads + (graph.get_num_masters()%nThreads==0? 0:1);
				const lvid_t begin = tid*chunkSize;
				const lvid_t end = (tid+1)*chunkSize < graph.get_num_masters()? (tid+1)*chunkSize : graph.get_num_masters();
				lvid_t vbuf[vbufSize];
				int nv=0;
				for(lvid_t lv=begin; lv<end; lv++){
					if(vertexCode.init(lv)){
						currSet.set(lv);
						vbuf[nv++]=lv;
						if(nv==vbufSize){
							currSet.insert(vbuf, nv);
							nv = 0;
						}
					}
				}
				currSet.insert(vbuf, nv);
				nv = 0;
			}
			//or
			if(syncState == PS_InitSet && tid==0){
				SourceVertexSet_t& srcSet = *(sourceVertexSet_ptr);
				const size_t n = srcSet.vertices.size();
				for(int i=0; i<n; i++){
					const vid_t v = srcSet.vertices[i];
					if(graph.get_owner(v)==procID){
						const lvid_t lv = graph.get_lvid_master(v);
						if(vertexCode.init(lv)){
							currSet.setAndInsert(lv);
						}
					}
				}
			}
			lc.barrier();//----------------------------------------3

			//Main-loop
			//for(int sstep = 0; sstep < max_nSupersteps; sstep++){
			while(1){
				lc.barrier();//----------------------------------------------4-0
				if(syncState == PS_Done) break;
				//1:Sync scatters
				//update native scatters
				///////////////////////////Phase-1////////////////////////////////
			//	assert(syncState == PS_Phase_Sync_Scatters);
				//1: send
				const int nThreads1=nThreads-1;//std::min(nThreads, nProcs);
				if(tid!=nThreads1)
				for(int p = tid; p < nProcs; p+=nThreads1){
					if(p==procID) continue;
					ScatterFormatedBuffer_t fbuf;
					byte* addr;
					do {
						addr = recvScatters(incomingQueue,currSet);
						if(!addr)
							addr = interBufPool.get_address();
					} while(!addr);//until get one. should not fail.
					fbuf.format(addr, interBufPool.bufSize());

					std::vector<lvid_t>& activeSet = currSet.exportSet(); 
					const int64_t beginIdx = 0;
					const int64_t endIdx = currSet.num_ready();
					if(dir==Out){
						Bitmap_t& scatterMap = graph.scatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(scatterMap.test(lv)){//lv has scatters in remote machine p				
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					} 
					if(dir==In){
						Bitmap_t& rscatterMap = graph.rscatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(rscatterMap.test(lv)){//lv has rscatters(combiners) in remote machine p				
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==InOut){
						Bitmap_t& scatterMap = graph.scatterMaps[p];
						Bitmap_t& rscatterMap = graph.rscatterMaps[p];
						for(int64_t idx = beginIdx; idx < endIdx; idx++){
							const lvid_t lv = activeSet[idx];
							if(scatterMap.test(lv)||rscatterMap.test(lv)){
								const vid_t v = graph.get_vid(lv);
								ScatterMsg_t msg;
								msg.dest = v;
								vertexCode.copy_scatter_data(msg.data, lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do{
										addr = recvScatters(incomingQueue,currSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}

					fbuf.getHeader()->setExtra(p);
					fbuf.getHeader()->setEnd();
					while(!outgoingQueue.enqueue(fbuf.address()))
						;//never fails! to ensure this.
				}
				//2:recv
				if(tid!=nThreads1)
				while(syncState == PS_Phase_Sync_Scatters){
					byte* addr = recvScatters(incomingQueue,currSet);
					if(addr){
						interBufPool.put_address(addr);
					} else 
						sched_yield();		
				}
				lc.barrier();//------------------------------------4-2

				//wait master to update currSet.
				//2:Compute
				///////////////////////////Phase-2: main Computation////////////////////////////////
				//assert(syncState == PS_Phase_Compute);
				if(tid==0){
					nFinishedGroups=0;
				}
				lc.barrier();//------------------------------------4-3
				do {
					//0
					for(int i=0; i<nGroups; i++){
						if(i!=ctx.groupID){
							byte* addr;
							do {
								addr = intraBufPool.get_address();
							} while(!addr);
							ctx.socketBuffers[i].format(addr, intraBufPool.bufSize());
						}
					}
					//1
					lvid_t buf[scheduleGrain];
					while(1){
						const int n = currSet.fetch(buf, scheduleGrain);
						if(n > 0){
							if(dir==Out || dir==InOut)
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								LocalGraph_t::edge_iterator it = localGraph.edge_begin(lu);
								const LocalGraph_t::edge_iterator end = localGraph.edge_end(lu);
								for(;it!=end; it++){
									const lvid_t lv = *it;
									vertexCode.scatter(ctx, lu, lv, it);
								}
							}
							if(dir==In || dir==InOut)
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								LocalGraph_t::redge_iterator it = localGraph.redge_begin(lu);
								const LocalGraph_t::redge_iterator end = localGraph.redge_end(lu);
								for(;it!=end; it++){
									const lvid_t lv = *it;
									vertexCode.scatter(ctx, lu, lv, it);
								}
							}
							for(int i=0; i<n; i++){
								const lvid_t lu = buf[i];
								if(vertexCode.assert_to_halt(lu))
									currSet.unset(lu);
							}
						} else
							break;
					}
					for(int i=0; i<nGroups; i++){
						if(i==ctx.groupID) continue;
						if(ctx.socketBuffers[i].count()>0){
							while(!group.putToChannel(i, ctx.socketBuffers[i].address()))
								;//wait and repeat
						} else {
							intraBufPool.put_address(ctx.socketBuffers[i].address());
						}
					}
					//2
					group.barrier();
					if(group.isLeader(tid)){
						AtomicFetchAdd(&nFinishedGroups, 1);
					}
					while(nFinishedGroups < nGroups){
						byte* addr;
						addr = recvFromChannel(ctx);
						if(addr){
							intraBufPool.put_address(addr);
						} else
							sched_yield();
					}
					while(1){
						byte* addr;
						addr = recvFromChannel(ctx);
						if(addr){
							intraBufPool.put_address(addr);
						} else
							break;
					}
				} while(0);
				lc.barrier();//------------------------------------4-4
				//3:Sync combiners
				///////////////////////////Phase-3: sync combiners////////////////////////////////
				//assert(syncState == PS_Phase_Sync_Combiners);
				if(tid!=nThreads1)
				for(int p = tid; p < nProcs; p+=nThreads1){
					if(p==procID) continue;
					CombineFormatedBuffer_t fbuf;
					byte* addr;
					do {
						addr = recvCombiners(ctx, incomingQueue, nextSet);
						if(!addr)
							addr = interBufPool.get_address();
					} while(!addr);//until get one. should not fail.
					fbuf.format(addr, interBufPool.bufSize());

					if(dir==Out){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& combinerMap = graph.combinerMaps[p];
						Bitmap_t::wordIterator it1 = combinerMap.begin();
						//const Bitmap_t::wordIterator end = activeBitset.getWordIterator(graph.get_end_lvid_combiner());
						//const Bitmap_t::wordIterator end1 = combinerMap.end();
						const lvid_t maxID = graph.get_end_lvid_combiner();
						for(lvid_t base = graph.get_begin_lvid_combiner();base < maxID; it++, it1++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==In){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& rcombinerMap = graph.rcombinerMaps[p];
						Bitmap_t::wordIterator it1 = rcombinerMap.begin();
						const lvid_t maxID = graph.get_end_lvid_scatter();
						for(lvid_t base = graph.get_begin_lvid_scatter();base < maxID; it++, it1++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}
					}
					if(dir==InOut){
						Bitmap_t& activeBitset = nextSet.exportBitset(); 
						Bitmap_t::wordIterator it = activeBitset.getWordIterator(graph.get_begin_lvid_combiner());
						Bitmap_t& combinerMap = graph.combinerMaps[p];
						Bitmap_t::wordIterator it1 = combinerMap.begin();
						Bitmap_t& rcombinerMap = graph.rcombinerMaps[p];
						Bitmap_t::wordIterator it2 = rcombinerMap.begin();
						const lvid_t maxID = graph.get_end_lvid_combiner();
						for(lvid_t base = graph.get_begin_lvid_combiner();base < maxID; it++, it1++, it2++, base+=64){//for
							lvid_t vbuf[64];
							const int n = util::interpretWord(vbuf, (*it)&(*it1 | *it2), base);
							for(int i=0; i<n; i++){
								const lvid_t lv = vbuf[i];
								const vid_t v = graph.get_vid(lv);
								CombineMsg_t msg;
								msg.dest = v;
								vertexCode.copy_combine_data(msg.data, lv);
								vertexCode.reset_combiner(lv);
								if(!fbuf.push_back(msg)){//full
									fbuf.getHeader()->setExtra(p);
									while(!outgoingQueue.enqueue(fbuf.address()))
										;//never fails! to ensure this.
									do {
										addr = recvCombiners(ctx, incomingQueue,nextSet);
										if(!addr)
											addr = interBufPool.get_address();
									} while(!addr);
									fbuf.format(addr, interBufPool.bufSize());
									fbuf.push_back(msg);									
								}
							}
						}

					}

					fbuf.getHeader()->setExtra(p);
					fbuf.getHeader()->setEnd();
					while(!outgoingQueue.enqueue(fbuf.address()))
						;//never fails! to ensure this.
				}
				//2:recv
				if(tid!=nThreads1)
				while(syncState == PS_Phase_Sync_Combiners){
					byte* addr = recvCombiners(ctx, incomingQueue,nextSet);
					if(addr){
						interBufPool.put_address(addr);
					} else 
						sched_yield();				
				}
				lc.barrier();//------------------------------------------4-5
				///////////////////////////Phase-4////////////////////////////////
				lc.barrier();//------------------------------------------4-6
				assert(syncState == PS_Phase_Apply);
				const size_t masterIdx_begin = graph.get_begin_lvid_master();
				const size_t masterIdx_end = graph.get_end_lvid_master();
				Bitmap_t& activeBitset = nextSet.exportBitset(); 
				const size_t nWords = util::align(masterIdx_end - masterIdx_begin)/64;
				//apply
				for(size_t i = masterIdx_begin/64+tid; i<nWords; i+=nThreads){
					uint64_t tmp = activeBitset.getWord(i);
					lvid_t lv = i*64;
					while(tmp){
						if(tmp & 0x1UL){
							if(vertexCode.apply(lv)) currSet.set(lv);
						}
						lv++;
						tmp>>=1;
					}
				}
				//refresh currSet for next superstep
				lvid_t vbuf[vbufSize];
				int nv = 0;
				Bitmap_t& activeBitset1 = currSet.exportBitset(); 
				for(size_t i = masterIdx_begin/64+tid; i<nWords; i+=nThreads){
					uint64_t tmp = activeBitset1.getWord(i);
					lvid_t lv = i*64;
					while(tmp){
						if(tmp & 0x1UL){
							vbuf[nv++] = lv;
							if(nv == vbufSize){
								currSet.insert(vbuf, nv);
								nv = 0;
							}
						}
						lv++;
						tmp>>=1;
					}
				}
				if(nv > 0){
					currSet.insert(vbuf, nv);
					nv = 0;
				}
				lc.barrier();//------------------------------------------4-7
			}
			return GRE_SUCCESS;
		}
	private:
		DistributedGraph_t& graph;
		SourceVertexSet_t* sourceVertexSet_ptr;
		//buffer pools
		parBufferPool_t intraBufPool;
		parBufferPool_t interBufPool;
		//active vertex set
		VertexSet_t currSet;
		VertexSet_t nextSet;
		//master-->............-->workers
		FFQueue_t incomingQueue;
		//workers-->............-->master
		FFQueue_t outgoingQueue;
		//Instance of Specific algorithmic primitives
		VertexCode& vertexCode;
		//
		bool isInitialized;
		//
		int superstep;
		int max_nSupersteps;
		//Thread State
		volatile CtrlState_t syncState;
		//Distributed control
		BufComm_t dc;
		//used only in thread group mode
		int nFinishedGroups;
		//timers
		Timer_t initVertexTimer; 
		Timer_t syncScatterTimer; 
		Timer_t syncCombinerTimer; 
		Timer_t localComputeTimer; 
		Timer_t localApplyTimer; 
		//Profiling
	public:
		double getInitVertexTime(){
			return initVertexTimer.getTime();
		}
		double getSyncScatterTime(){
			return syncScatterTimer.getTime();
		}
		double getSyncCombinerTime(){
			return syncCombinerTimer.getTime();
		}
		double getLocalComputeTime(){
			return localComputeTimer.getTime();
		}
		double getLocalApplyTime(){
			return localApplyTimer.getTime();
		}
		int getSuperstepCount(){
			return superstep;
		}
	private:
		//
		//handle received combiners
		template <typename Ctx_t>
		inline byte* recvCombiners(Ctx_t& ctx, FFQueue_t& inComingQueue, VertexSet_t& vSet){
			byte* addr = inComingQueue.dequeue();
			if(addr){
				CombineFormatedBuffer_t buf(addr, interBufPool.bufSize());
				typename CombineFormatedBuffer_t::data_iterator it = buf.begin();
				const typename CombineFormatedBuffer_t::data_iterator end = buf.end();
				for(; it!=end; it++){
					const lvid_t lv = graph.get_lvid_master(it->dest);
					//sendMessage(ctx, lv, it->data);
					typename Ctx_t::Lock_t& lock = ctx.vlock.getLockRef(lv);
					if(vertexCode.combine(lv, it->data, lock)){
						vSet.testAndSet(lv);
					}
				}
			} 
			return addr;
		}
		//
		//handle received scatters
		inline byte* recvScatters(FFQueue_t& inComingQueue, VertexSet_t& vSet){
			byte* addr = inComingQueue.dequeue();
			if(addr){
				ScatterFormatedBuffer_t buf(addr, interBufPool.bufSize());
				lvid_t vbuf[vbufSize];
				int n=0;
				typename ScatterFormatedBuffer_t::data_iterator it = buf.begin();
				const typename ScatterFormatedBuffer_t::data_iterator end = buf.end();
				for(; it!=end; it++){
					const lvid_t lv = graph.get_lvid_scatter(it->dest);
					vertexCode.update_scatter(lv, it->data);
					//if(!vSet.testAndSet(lv)){//always succeed!
						vbuf[n++]=lv;
						if(n==vbufSize){//full
							vSet.insert(vbuf, vbufSize);
							n=0;
						}
					//}
				}
				if(n>0)
					vSet.insert(vbuf, n);
			} 
			return addr;
		}
		//
		//
		inline byte* recvFromChannel(TG_Context_t& ctx){
			byte* addr = ctx.group.getFromChannel();
			if(addr){
				SocketFormatedBuffer_t buf(addr, intraBufPool.bufSize());
				typename SocketFormatedBuffer_t::data_iterator it = buf.begin();
				const typename SocketFormatedBuffer_t::data_iterator end = buf.end();
				for(; it!=end; it++){
					const lvid_t lv = it->dest;
					typename TG_Context_t::Lock_t& lock = ctx.g_vlock.getLockRef(lv);
					if(vertexCode.combine(lv, it->data, lock)){
						ctx.nextSet.testAndSet(lv);
					}
				}
			} 
			return addr;
		}
	public:
		//send message
		template <typename Ctx_t, typename Data_t>
		inline void sendMessage(Ctx_t& ctx, lvid_t lv, Data_t& data){
			typename Ctx_t::Lock_t& lock = ctx.vlock.getLockRef(lv);
			if(vertexCode.combine(lv, data, lock))
				ctx.nextSet.testAndSet(lv);
		}

		template <typename Data_t>
		inline void sendMessage(TP_Context_t& ctx, lvid_t lv, Data_t& data){
			typename TP_Context_t::Lock_t& lock = ctx.vlock.getLockRef(lv);
			if(vertexCode.combine(lv, data, lock))
				ctx.nextSet.testAndSet(lv);
		}
		template <typename Data_t>
		inline void sendMessage(TG_Context_t& ctx, lvid_t lv, Data_t& data){
			const int vgid = lv % ctx.nGroups;
			if(vgid == ctx.groupID){
				typename TG_Context_t::Lock_t& lock = ctx.g_vlock.getLockRef(lv);//group vlock
				if(vertexCode.combine(lv, data, lock))
					ctx.nextSet.testAndSet(lv);
			} else {
				SocketMsg_t msg(lv, data);
				if(!ctx.socketBuffers[vgid].push_back(msg)){
					while(!ctx.group.putToChannel(vgid, ctx.socketBuffers[vgid].address()))
						;//
					byte* addr;
					do {
						addr = recvFromChannel(ctx);
						if(!addr)
							addr = intraBufPool.get_address();
					} while(!addr);
					ctx.socketBuffers[vgid].format(addr, intraBufPool.bufSize());
					ctx.socketBuffers[vgid].push_back(msg);
				}
			}
		}
	private:
		template<typename DC_T>
		void output_statistics(DC_T& dc){
			vertexCode.output_statistics(dc);
		}
	}; //end SyncEngine
}//end namespace
#endif
