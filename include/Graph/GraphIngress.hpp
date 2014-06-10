#ifndef GRE_Graph_GraphIngress_HPP
#define GRE_Graph_GraphIngress_HPP

#include "Types.hpp"
#include "Graph/LocalGraph.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Util/Util.hpp"
#include "Util/BufferPool.hpp"
#include "Util/Buffer.hpp"
#include "Util/Bitmap.hpp"
#include "Util/HashIndex.hpp"
#include "Comm/BufferingComm.hpp"
#include "Comm/Protocol.hpp"
#include "Comm/CommUnderlying.hpp"
#include "Config/UsrConfig.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <cstdio>
namespace GRE{
	namespace Graph{
		//Misc
		struct IngressParameters{
			//internal fields
			std::string format;
			std::string method;
			std::string prefix;
			//routines
			IngressParameters():prefix(),format("ascii"),method("edgelists"){}
			IngressParameters(std::string sprefix):prefix(sprefix),format("ascii"),method("edgelists"){}
			IngressParameters(std::string sprefix, std::string sformat, std::string smethod):prefix(sprefix),format(sformat),method(smethod){}
			void setFormat(std::string sformat){
				if(sformat=="bin") format=sformat;
				else format="ascii";
			}
			void setMethod(std::string& smethod){
				if(smethod=="partitions") method=smethod;
				else method="edgelists";
			}
			void setPrefix(std::string& sprefix){
				prefix=sprefix;
			}
			std::string& getPrefix()  {return prefix;}
			bool isBinary() const {return format=="bin";}
			bool isAscii() const {return format=="ascii";}
			bool isEdgelists() const {return method=="edgelists";}
			bool isPartitions() const {return method=="partitions";}
		};
		//Graph Ingress Infras
		template<typename DistributedGraph_t>
		class GraphIngress{
		public:
			typedef typename GRE::COMM::Underlying Underlying_t;
			typedef typename GRE::COMM::Protocol::Header Header_t;

			typedef typename GRE::Util::BufferPool BufferPool_t;
			typedef typename GRE::Util::FormatedBuffer<Header_t, edge_t> FormatedBuffer_t;
			typedef typename GRE::Util::FormatedBuffer<Header_t, vid_t> FormatedBuffer1_t;
			typedef typename GRE::Util::Bitmap Bitmap_t;
			typedef typename GRE::Util::FixedBitmap<GRE::MAX_NUM_MACHINES> FixedBitmap_t;
			typedef typename GRE::Util::HashIndex<vid_t, lvid_t> HashIndex_t;

			typedef typename GRE::COMM::BufComm<BufferPool_t> BufComm_t;
			typedef typename DistributedGraph_t::LocalGraph_t LocalGraph_t;
		public:
			//
			GraphIngress(DistributedGraph_t& dg_r,Underlying_t& uc_r, BufferPool_t& bp_r)
				:distributedGraph(dg_r), comm(uc_r), bufComm(uc_r, bp_r), bufPool(bp_r) {}
			~GraphIngress(){
				//bufComm.finalize();
			}
			//
			//Graph Paritioning Inside
			//
			template<Dir dir>
			bool ingress(IngressParameters& params){
				//0: simply load meta data, like nV, nE
				loadMetaFromFile(params.prefix);
				const int nProcs = comm.size();
				const int procID = comm.rank();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nEdges = distributedGraph.gStat_nEdges;
				const size_t nMasters = nVertices/nProcs + (procID < (nVertices%nProcs)? 1:0);
				distributedGraph.gStat_nSubgraphs = nProcs;
				//1: load graph topology
				edges.reserve(1.05*nEdges/nProcs);//release factor: 0.05
				if(procID == 0) std::cout << "Use load method "<<params.method<<std::endl;
				if(params.method == "partition"){
					if(dir==Out)
						loadEdgePartitions(params.prefix);
					else {//not-support
						std::cerr<<"Error: don't support loading partitions of In edges."<<std::endl;
						return false;
					}
				} else { //partition in place
					//load from text and graph partitioning --> edge lists
					if(params.format == "bin")
						loadEdgesFromBinaryFiles(params.prefix);
					else {
						loadEdgesFromFiles(params.prefix);
						/*
						if(dir==InOut)
							loadEdgesFromFiles1(params.prefix);
						else
							loadEdgesFromFiles(params.prefix);
						*/
					}
				}

				//2: build VertexBitmap
				Bitmap_t scatterBitmap(nVertices);
				Bitmap_t combinerBitmap(nVertices);
				size_t nScatters = 0;
				size_t nCombiners = 0;
				size_t nBiAgents = 0;
				const size_t n = edges.size();
				#pragma omp parallel 
				{
					size_t myNScatters = 0;
					size_t myNCombiners = 0;
					size_t myNBiAgents = 0;
					#pragma omp for
					for(size_t i=0; i<n; i++)
					{
						const vid_t src = edges[i].first;
						const vid_t dst = edges[i].second;
						if(getPartID(src, nProcs)!=procID){
							if(!scatterBitmap.testAndSet(src))
								myNScatters++;
						}
						if(getPartID(dst, nProcs)!=procID){
							if(!combinerBitmap.testAndSet(dst))
								myNCombiners++;
						}
					}
					AtomicFetchAdd(&nScatters, myNScatters);
					AtomicFetchAdd(&nCombiners, myNCombiners);
					//
					#pragma omp for
					for(size_t i=0; i<nVertices;i++){
						if(scatterBitmap.test(i) && combinerBitmap.test(i))
							myNBiAgents++;
					}
					AtomicFetchAdd(&nBiAgents, myNBiAgents);
				}
				//std::cerr<<"nBiAgents: "<<nBiAgents<<std::endl;
				//std::cerr<<"nCombiners: "<<nCombiners<<std::endl;
				//std::cerr<<"nScatters: "<<nScatters<<std::endl;
				//std::cerr<<"nMasters: "<<nMasters<<std::endl;
				std::cerr<<procID<<"nEdges: "<<n<<std::endl;
				if(dir==In || dir==InOut){
					nCombiners = nCombiners+nScatters-nBiAgents;
					nScatters = nCombiners;
				}
				//3: build lvid2vid
				//M+C
				std::vector<vid_t>& lvid2vid = distributedGraph.lvid2vid;
				lvid2vid.reserve(GRE::Util::align(nMasters)+nCombiners);
				//add masters
				for(vid_t i=procID; i<nVertices; i+=nProcs){
					lvid2vid.push_back(i);		
				}
				assert(lvid2vid.size()==nMasters);
				//padding...
				for(vid_t i=nMasters; i<GRE::Util::align(nMasters); i++){
					lvid2vid.push_back(-1);//-1 is invalid vertex id.		
				}

				if(dir==Out)//add combiners
				for(vid_t i=0; i<nVertices; i++){
					if(combinerBitmap.test(i))
						lvid2vid.push_back(i);
				}
				else//add agents
				for(vid_t i=0; i<nVertices; i++){
					if(combinerBitmap.test(i) || scatterBitmap.test(i))
						lvid2vid.push_back(i);
				}
				assert(lvid2vid.size()==GRE::Util::align(nMasters)+nCombiners);

				//build combinermaps for each remote process.
				std::vector<Bitmap_t>& combinerMaps = distributedGraph.combinerMaps;
				combinerMaps.resize(nProcs);
				for(int i=0; i<nProcs; i++){
					if(i != procID){
						combinerMaps[i].resize(nCombiners);//Shoud give more guarrantte on bounding.
					}
				}
				do{
					const size_t base = GRE::Util::align(nMasters);
					for(size_t i=0; i<nCombiners; i++){
						const vid_t v = lvid2vid[base+i];
						const int owner = getPartID(v, nProcs);
						combinerMaps[owner].set_serial(i);
					}
				} while(0);

				//build rcombinermaps for each remote process.
				if(dir!=Out){
					std::vector<Bitmap_t>& rcombinerMaps = distributedGraph.rcombinerMaps;
					rcombinerMaps.resize(nProcs);
					for(int i=0; i<nProcs; i++){
						if(i != procID){
							rcombinerMaps[i].resize(nScatters);//Shoud give more guarrantte on bouding.
						}
					}
					do{
						const size_t base = GRE::Util::align(nMasters);
						for(size_t i=0; i<nScatters; i++){
							const vid_t v = lvid2vid[base+i];
							const int owner = getPartID(v, nProcs);
							rcombinerMaps[owner].set_serial(i);
						}
					} while(0);
				}

				//4: build vid2lvid
				//M+S
				HashIndex_t& vtxIndex = distributedGraph.vtxIndex;	
				if(dir==Out){
					std::vector<vid_t> lvid2vid1;
					lvid2vid1.reserve(GRE::Util::align(nMasters)+nScatters);
					//add masters
					for(vid_t i=procID; i<nVertices; i+=nProcs){
						lvid2vid1.push_back(i);		
					}
					assert(lvid2vid1.size()==nMasters);
					//padding...
					for(vid_t i=nMasters; i<GRE::Util::align(nMasters); i++){
						lvid2vid1.push_back(-1);//-1 is invalid vertex id.		
					}
					//add scatters
					for(vid_t i=0; i<nVertices; i++){
						if(scatterBitmap.test(i))
							lvid2vid1.push_back(i);
					}
					assert(lvid2vid1.size()==GRE::Util::align(nMasters)+nScatters);
					if(nScatters>0)
						vtxIndex.buildFromArray(lvid2vid1, GRE::Util::align(nMasters), nScatters);
					else std::cerr<<procID<<"hash index build....nscatters is "<<nScatters<<std::endl;
					lvid2vid1.clear();
				} else {
					if(nScatters>0)
						vtxIndex.buildFromArray(lvid2vid, GRE::Util::align(nMasters), nScatters);
					else std::cerr<<procID<<"hash index build....nscatters is "<<nScatters<<std::endl;
				}

				//5: build CSR
				//
				std::vector<ledge_t> ledges;
				ledges.resize(n);

				if(dir==Out){
					HashIndex_t combinerIndex(-1);	
					if(nCombiners > 0)
						combinerIndex.buildFromArray(lvid2vid, GRE::Util::align(nMasters), nCombiners);
				
					#pragma omp parallel for
					for(size_t i=0; i<n; i++) {
						const vid_t src = edges[i].first;
						const vid_t dst = edges[i].second;
						if(src%nProcs == procID)
							ledges[i].first = src/nProcs;
						else
							ledges[i].first = vtxIndex.value(src);

						if(dst%nProcs == procID)
							ledges[i].second = dst/nProcs;
						else
							ledges[i].second = combinerIndex.value(dst);
					}
					combinerIndex.clear();
				} else {
					//
					#pragma omp parallel for
					for(size_t i=0; i<n; i++) {
						const vid_t src = edges[i].first;
						const vid_t dst = edges[i].second;
						if(src%nProcs == procID)
							ledges[i].first = src/nProcs;
						else
							ledges[i].first = vtxIndex.value(src);

						if(dst%nProcs == procID)
							ledges[i].second = dst/nProcs;
						else
							ledges[i].second = vtxIndex.value(dst);
					}
				}
				
				distributedGraph.graph.build(ledges, GRE::Util::align(nMasters)+nScatters);
				if(dir==In || dir==InOut)
					distributedGraph.graph.rbuild(ledges, GRE::Util::align(nMasters)+nCombiners);
				ledges.clear();
				//6:
				distributedGraph.lStat_nMasters = nMasters;
				distributedGraph.lStat_nScatters = nScatters;
				distributedGraph.lStat_nCombiners = nCombiners;
				distributedGraph.lStat_nAgents = nCombiners+nScatters-nBiAgents;
				//7:
				if(dir!=Out)
					syncrScatters();
			}

			bool ingress(const IngressParameters& params){
				ingress<Out>(params.prefix, params.format, params.method);
			}
			bool finalize(){
				//ok, I am ready!
				distributedGraph.isReady = true;
				edges.clear();
			}
		private:
			//Infra
			DistributedGraph_t& distributedGraph;
			Underlying_t& comm;
			BufComm_t bufComm;
			BufferPool_t& bufPool;
			//edge lists in memory during loading
			std::vector<edge_t> edges;
		private:
			void loadMetaFromFile(const std::string& prefix){
				std::string filename = prefix + ".meta";
				std::ifstream infile(filename.c_str());
				infile >> distributedGraph.gStat_nVertices;
				infile >> distributedGraph.gStat_nEdges;
			}
			void loadEdgesFromBinaryFiles(const std::string& prefix){
				//open stream
				//read-in
				const int partID = comm.rank();
				const int nParts = comm.size();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nEdges = distributedGraph.gStat_nEdges;
				const size_t nMasters = nVertices/nParts + (nVertices%nParts > partID? 1:0);

				std::vector<FormatedBuffer_t> buffers;
				buffers.resize(nParts);
				for(int i=0;i<nParts;i++){
					if(i==partID) continue;
					const byte* addr = bufPool.get_address();//assert(addr);
					const size_t capacity = bufPool.bufSize();
					buffers[i].format(addr, capacity);
				}

				std::vector<FixedBitmap_t> srcMap;
				std::vector<FixedBitmap_t> dstMap;
				std::vector<size_t> edgeCounts;

				srcMap.resize(nVertices);
				dstMap.resize(nVertices);
				edgeCounts.resize(nParts, 0);
			
				//TODO:Open or close the pre-decision.
				/*
				#pragma omp parallel for
				for(vid_t i=0; i<nVertices; i++){
					srcMap[i].set_serial(getPartID(i, nParts));
				}
				#pragma omp parallel for
				for(vid_t i=0; i<nVertices; i++){
					dstMap[i].set_serial(getPartID(i, nParts));
				}
				*/

				comm.barrier();
				
				if(!bufComm.launchAllRecv())
					std::cerr << "Fatal Err: fail to launch recv." <<std::endl;

				EdgeDecision edgeDecision(nParts, partID, srcMap, dstMap, edgeCounts);
				std::vector<Bitmap_t>& scatterMaps = distributedGraph.scatterMaps;
				std::vector<int>& outdegrees = distributedGraph.outdegrees;
				scatterMaps.resize(nParts);
				outdegrees.resize(nMasters,0);
				for(int i=0; i<nParts; i++){
					if(i != partID){
						scatterMaps[i].resize(nMasters);//Shoud give more guarrantte on bouding.
					}
				}
				for(int i=partID; ;i+=nParts)
				{
					std::string filename(prefix);
					filename = filename + ".edges."+ GRE::Util::toString(i);
					std::ifstream infile(filename.c_str(),std::ios_base::in | std::ios_base::binary);
					if(!infile) {
						std::cerr<<partID<<" : Ok! Finish loading edgelist files."<<std::endl;
						break;
					}
					std::cout<<"Loading binary edges from"<<filename<<"..."<<std::endl;

					while(infile.good() && !infile.eof()){
						vid_t src, dest;
						if(!(infile >> src) || !(infile >> dest)){
							std::cerr<<"Fetal Err in reading edge list files"<<std::endl;
							exit(-1);
						}
						int epid = edgeDecision.getPartID(src, dest);			
						const edge_t edge(src, dest);
						outdegrees[src/nParts]++;
						if(epid == partID){
							edges.push_back(edge);
						} else {
							//compute scatterMap
							scatterMaps[epid].set_serial(src/nParts);
							//send in buffering mode
							if(!buffers[epid].push_back(edge)){
								buffers[epid].getHeader()->setExtra(partID);//rank
								bufComm.isend(buffers[epid], epid);
								byte* addr;
								do {
									bool extra = false;
									bufComm.checkAllSend();
									bufComm.checkAllRecv(extra);
									//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-1."<<std::endl;
									typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
									byte* addr1;
									while(queue.get(addr1)){
										FormatedBuffer_t buf;
										buf.open(addr1, bufPool.bufSize());
										if(buf.getHeader()->isEnd()){
											bufComm.cancelRecv(buf.getHeader()->getExtra());
										}
										FormatedBuffer_t::data_iterator cur = buf.begin();
										const FormatedBuffer_t::data_iterator end = buf.end();
										for(; cur != end; cur++){
											edges.push_back(*cur);
										}
										buf.close();
										bufPool.put_address(addr1);
										if(!extra) {
											bufComm.resumeRecv();
											extra = true;
										}
									}
									addr = bufPool.get_address();
								} while(!addr);
								buffers[epid].format(addr,bufPool.bufSize());
								buffers[epid].push_back(edge);	
							}
						}
					}//end-while-reading
					do {
						bool extra = false;
						bufComm.checkAllSend();
						bufComm.checkAllRecv(extra);
						//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-2."<<std::endl;
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();	
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					} while(0);
				}//end-for-reading-all
				//Broadcast to others: I am finished.
				for(int i=0; i<nParts; i++){
					if(i!=partID){
						buffers[i].getHeader()->setEnd();
						buffers[i].getHeader()->setExtra(partID);
						bufComm.isend(buffers[i], i);
					}
				}
				do{
					bool extra = false;
					while(bufComm.checkAllRecv(extra)){
						bufComm.checkAllSend();
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
	
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					}
				} while(0);
				while(bufComm.checkAllSend());
				
				//optional sync
				comm.barrier();
				srcMap.clear();
				dstMap.clear();
				edgeCounts.clear();
			}
			
			void loadEdgesFromFiles(const std::string& prefix){
				//open stream
				//read-in
				const int partID = comm.rank();
				const int nParts = comm.size();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nEdges = distributedGraph.gStat_nEdges;
				const size_t nMasters = nVertices/nParts + (nVertices%nParts > partID? 1:0);

				std::vector<FormatedBuffer_t> buffers;
				buffers.resize(nParts);
				for(int i=0;i<nParts;i++){
					if(i==partID) continue;
					const byte* addr = bufPool.get_address();//assert(addr);
					const size_t capacity = bufPool.bufSize();
					buffers[i].format(addr, capacity);
				}

				std::vector<FixedBitmap_t> srcMap;
				std::vector<FixedBitmap_t> dstMap;
				std::vector<size_t> edgeCounts;

				srcMap.resize(nVertices);
				dstMap.resize(nVertices);
				edgeCounts.resize(nParts, 0);
		
				/*
				//#ifdef OPTIMISE_GP
				//#pragma omp parallel for
				for(vid_t i=0; i<1024; i++){
					int j=(rand()+partID)%nVertices;
					srcMap[j].set_serial(getPartID(j, nParts));
					//srcMap[i].set_serial(getPartID(i, nParts));
				}
				//#pragma omp parallel for
				for(vid_t i=0; i<1024; i++){
					int j=(rand()+partID)%nVertices;
					//dstMap[i].set_serial(getPartID(i, nParts));
					dstMap[j].set_serial(getPartID(j, nParts));
				}
				//#endif
				*/

				comm.barrier();
				
				if(!bufComm.launchAllRecv())
					std::cerr << "Fatal Err: fail to launch recv." <<std::endl;

				//GreedyEdgeDecision edgeDecision(nParts, partID, srcMap, dstMap, edgeCounts);
				LDGEdgeDecision edgeDecision(nParts, partID, nEdges/nParts/nParts, srcMap, dstMap, edgeCounts);
				std::vector<Bitmap_t>& scatterMaps = distributedGraph.scatterMaps;
				std::vector<int>& outdegrees = distributedGraph.outdegrees;
				scatterMaps.resize(nParts);
				outdegrees.resize(nMasters,0);
				for(int i=0; i<nParts; i++){
					if(i != partID){
						scatterMaps[i].resize(nMasters);//Shoud give more guarrantte on bouding.
					}
				}
				for(int i=partID; ;i+=nParts)
				{
					std::string filename(prefix);
					filename = filename + ".edges."+ GRE::Util::toString(i);
					std::ifstream infile(filename.c_str());
					if(!infile) {
						std::cerr<<partID<<" : Ok! Finish loading edgelist files."<<std::endl;
						break;
					}
					std::cout<<"Loading edges from"<<filename<<"...";
					//SNAP and TSV
					while(infile.good() && !infile.eof()){
						std::string line;
						std::getline(infile, line);
						if(line.empty() || line[0]=='#') continue;
						if(infile.fail()) break;
						char* tmpptr;
						const vid_t src = strtoul(line.c_str(), &tmpptr, 10);
						const vid_t dest = strtoul(tmpptr, NULL, 10);
						int epid = edgeDecision.getPartID(src, dest);			
						const edge_t edge(src, dest);
						outdegrees[src/nParts]++;
						if(epid == partID){
							edges.push_back(edge);
						} else {
							//compute scatterMap
							scatterMaps[epid].set_serial(src/nParts);
							//send in buffering mode
							if(!buffers[epid].push_back(edge)){
								buffers[epid].getHeader()->setExtra(partID);//rank
								bufComm.isend(buffers[epid], epid);
								byte* addr;
								do {
									bool extra = false;
									bufComm.checkAllSend();
									bufComm.checkAllRecv(extra);
									//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-1."<<std::endl;
									typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
									byte* addr1;
									while(queue.get(addr1)){
										FormatedBuffer_t buf;
										buf.open(addr1, bufPool.bufSize());
										if(buf.getHeader()->isEnd()){
											bufComm.cancelRecv(buf.getHeader()->getExtra());
										}
										FormatedBuffer_t::data_iterator cur = buf.begin();
										const FormatedBuffer_t::data_iterator end = buf.end();
										for(; cur != end; cur++){
											edges.push_back(*cur);
										}
										buf.close();
										bufPool.put_address(addr1);
										if(!extra) {
											bufComm.resumeRecv();
											extra = true;
										}
									}
									addr = bufPool.get_address();
								} while(!addr);
								buffers[epid].format(addr,bufPool.bufSize());
								buffers[epid].push_back(edge);	
							}
						}
					}//end-while-reading
					do {
						bool extra = false;
						bufComm.checkAllSend();
						bufComm.checkAllRecv(extra);
						//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-2."<<std::endl;
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();	
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					} while(0);
					std::cout<<"ok!"<<std::endl;
				}//end-for-reading-all
				//Broadcast to others: I am finished.
				for(int i=0; i<nParts; i++){
					if(i!=partID){
						buffers[i].getHeader()->setEnd();
						buffers[i].getHeader()->setExtra(partID);
						bufComm.isend(buffers[i], i);
					}
				}
				/*
				do{
					bool extra = false;
					while(bufComm.checkAllRecv(extra)){
						bufComm.checkAllSend();
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
	
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					}
				} while(0);
				*/
				do{
					while(bufComm.checkAllRecv()){
						bufComm.checkAllSend();
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
	
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
						}
					}
				} while(0);

				while(bufComm.checkAllSend());
				
				//optional sync
				comm.barrier();
				srcMap.clear();
				dstMap.clear();
				edgeCounts.clear();
			}
			void loadEdgesFromFiles1(const std::string& prefix){
				//open stream
				//read-in
				const int partID = comm.rank();
				const int nParts = comm.size();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nEdges = distributedGraph.gStat_nEdges;
				const size_t nMasters = nVertices/nParts + (nVertices%nParts > partID? 1:0);

				std::vector<FormatedBuffer_t> buffers;
				buffers.resize(nParts);
				for(int i=0;i<nParts;i++){
					if(i==partID) continue;
					const byte* addr = bufPool.get_address();//assert(addr);
					const size_t capacity = bufPool.bufSize();
					buffers[i].format(addr, capacity);
				}

				std::vector<FixedBitmap_t> vtxMap;
				std::vector<size_t> edgeCounts;

				vtxMap.resize(nVertices);
				edgeCounts.resize(nParts, 0);
			
				//TODO:Open or close the pre-decision.
				//#pragma omp parallel for
				//for(vid_t i=0; i<nVertices; i++){
				//#ifdef OPTIMISE_GP
				for(vid_t i=0; i<2048; i++){
					int j=(rand()+partID)%nVertices;
					vtxMap[j].set_serial(getPartID(j, nParts));
				}
				//#endif
				
				comm.barrier();
				
				if(!bufComm.launchAllRecv())
					std::cerr << "Fatal Err: fail to launch recv." <<std::endl;

				EdgeDecision1 edgeDecision(nParts, partID, vtxMap, edgeCounts);
				std::vector<Bitmap_t>& scatterMaps = distributedGraph.scatterMaps;
				std::vector<int>& outdegrees = distributedGraph.outdegrees;
				scatterMaps.resize(nParts);
				outdegrees.resize(nMasters,0);
				for(int i=0; i<nParts; i++){
					if(i != partID){
						scatterMaps[i].resize(nMasters);//Shoud give more guarrantte on bouding.
					}
				}
				for(int i=partID; ;i+=nParts)
				{
					std::string filename(prefix);
					filename = filename + ".edges."+ GRE::Util::toString(i);
					std::ifstream infile(filename.c_str());
					if(!infile) {
						std::cerr<<partID<<" : Ok! Finish loading edgelist files."<<std::endl;
						break;
					}
					std::cout<<"Loading edges from"<<filename<<"..."<<std::endl;
					//SNAP and TSV
					while(infile.good() && !infile.eof()){
						std::string line;
						std::getline(infile, line);
						if(line.empty() || line[0]=='#') continue;
						if(infile.fail()) break;
						char* tmpptr;
						const vid_t src = strtoul(line.c_str(), &tmpptr, 10);
						const vid_t dest = strtoul(tmpptr, NULL, 10);
						int epid = edgeDecision.getPartID(src, dest);			
						const edge_t edge(src, dest);
						outdegrees[src/nParts]++;
						if(epid == partID){
							edges.push_back(edge);
						} else {
							//compute scatterMap
							scatterMaps[epid].set_serial(src/nParts);
							//send in buffering mode
							if(!buffers[epid].push_back(edge)){
								buffers[epid].getHeader()->setExtra(partID);//rank
								bufComm.isend(buffers[epid], epid);
								byte* addr;
								do {
									bool extra = false;
									bufComm.checkAllSend();
									bufComm.checkAllRecv(extra);
									//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-1."<<std::endl;
									typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
									byte* addr1;
									while(queue.get(addr1)){
										FormatedBuffer_t buf;
										buf.open(addr1, bufPool.bufSize());
										if(buf.getHeader()->isEnd()){
											bufComm.cancelRecv(buf.getHeader()->getExtra());
										}
										FormatedBuffer_t::data_iterator cur = buf.begin();
										const FormatedBuffer_t::data_iterator end = buf.end();
										for(; cur != end; cur++){
											edges.push_back(*cur);
										}
										buf.close();
										bufPool.put_address(addr1);
										if(!extra) {
											bufComm.resumeRecv();
											extra = true;
										}
									}
									addr = bufPool.get_address();
								} while(!addr);
								buffers[epid].format(addr,bufPool.bufSize());
								buffers[epid].push_back(edge);	
							}
						}
					}//end-while-reading
					do {
						bool extra = false;
						bufComm.checkAllSend();
						bufComm.checkAllRecv(extra);
						//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-2."<<std::endl;
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();	
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					} while(0);
				}//end-for-reading-all
				//Broadcast to others: I am finished.
				for(int i=0; i<nParts; i++){
					if(i!=partID){
						buffers[i].getHeader()->setEnd();
						buffers[i].getHeader()->setExtra(partID);
						bufComm.isend(buffers[i], i);
					}
				}
				/*
				do{
					bool extra = false;
					while(bufComm.checkAllRecv(extra)){
						bufComm.checkAllSend();
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
	
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
							if(!extra) {
								bufComm.resumeRecv();
								extra = true;
							}
						}
					}
				} while(0);
				*/
				do{
					while(bufComm.checkAllRecv()){
						bufComm.checkAllSend();
						typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
						byte* addr;
						while(queue.get(addr)){
							FormatedBuffer_t buf(addr, bufPool.bufSize());
							if(buf.getHeader()->isEnd()){
								bufComm.cancelRecv(buf.getHeader()->getExtra());
							}
	
							FormatedBuffer_t::data_iterator cur = buf.begin();
							const FormatedBuffer_t::data_iterator end = buf.end();
							for(; cur != end; cur++){
								edges.push_back(*cur);
							}
							bufPool.put_address(addr);
						}
					}
				} while(0);

				while(bufComm.checkAllSend());
				
				//optional sync
				comm.barrier();
				vtxMap.clear();
				edgeCounts.clear();
			}
			void syncrScatters(){
				const int partID = comm.rank();
				const int nParts = comm.size();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nMasters = nVertices/nParts + (nVertices%nParts > partID? 1:0);
				const size_t nCombiners = distributedGraph.lStat_nCombiners;
				
				std::vector<FormatedBuffer1_t> buffers;
				buffers.resize(nParts);
				for(int i=0;i<nParts;i++){
					if(i==partID) continue;
					const byte* addr = bufPool.get_address();//assert(addr);
					const size_t capacity = bufPool.bufSize();
					buffers[i].format(addr, capacity);
				}
				if(!bufComm.launchAllRecv())
					std::cerr << "Fatal Err: fail to launch recv." <<std::endl;


				std::vector<Bitmap_t>& rscatterMaps = distributedGraph.rscatterMaps;
				rscatterMaps.resize(nParts);
				for(int i=0; i<nParts; i++){
					if(i != partID){
						rscatterMaps[i].resize(nMasters);//Shoud give more guarrantte on bouding.
					}
				}
				const size_t base = GRE::Util::align(nMasters);
				std::vector<vid_t>& lvid2vid = distributedGraph.lvid2vid;
				for(size_t i=0; i<nCombiners; i++){
					const vid_t v = lvid2vid[base+i];
					const int owner = getPartID(v, nParts);
					if(!buffers[owner].push_back(v)){
						buffers[owner].getHeader()->setExtra(partID);//rank
						bufComm.isend(buffers[owner], owner);
						byte* addr;
						do {
							bool extra = false;
							bufComm.checkAllSend();
							bufComm.checkAllRecv(extra);
							//if(!bufComm.checkAllRecv()) std::cerr<<partID<<" has finished-1."<<std::endl;
							typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
							byte* addr1;
							while(queue.get(addr1)){
								FormatedBuffer1_t buf;
								buf.open(addr1, bufPool.bufSize());
								if(buf.getHeader()->isEnd()){
									bufComm.cancelRecv(buf.getHeader()->getExtra());
								}
								const int pid=buf.getHeader()->getExtra();
								FormatedBuffer1_t::data_iterator cur = buf.begin();
								const FormatedBuffer1_t::data_iterator end = buf.end();
								for(; cur != end; cur++){
									rscatterMaps[pid].set_serial((*cur)/nParts);
								}
								buf.close();
								bufPool.put_address(addr1);
								if(!extra) {
									bufComm.resumeRecv();
									extra = true;
								}
							}
							addr = bufPool.get_address();
						} while(!addr);
						buffers[owner].format(addr,bufPool.bufSize());
						buffers[owner].push_back(v);	
					}
				}
				for(int i=0; i<nParts; i++){
					if(i!=partID){
						buffers[i].getHeader()->setEnd();
						buffers[i].getHeader()->setExtra(partID);
						bufComm.isend(buffers[i], i);
					}
				}

				while(bufComm.checkAllRecv()){
					bufComm.checkAllSend();
					typename BufComm_t::incomingQueue_t& queue = bufComm.getIncomingQueue();
					byte* addr;
					while(queue.get(addr)){
						FormatedBuffer1_t buf(addr, bufPool.bufSize());
						if(buf.getHeader()->isEnd()){
							bufComm.cancelRecv(buf.getHeader()->getExtra());
						}
						const int pid=buf.getHeader()->getExtra();
						FormatedBuffer1_t::data_iterator cur = buf.begin();
						const FormatedBuffer1_t::data_iterator end = buf.end();
						for(; cur != end; cur++){
							rscatterMaps[pid].set_serial((*cur)/nParts);
						}
						bufPool.put_address(addr);
					}
				}

				while(bufComm.checkAllSend());
				
				//optional sync
				comm.barrier();
			}
			void loadEdgePartitions(const std::string& prefix){
				const int partID = comm.rank();
				const int nParts = comm.size();
				const size_t nVertices = distributedGraph.gStat_nVertices;
				const size_t nEdges = distributedGraph.gStat_nEdges;
				const size_t nMasters = nVertices/nParts + (nVertices%nParts > partID? 1:0);
				//1
				do {
					std::vector<int>& outdegrees = distributedGraph.outdegrees;
					outdegrees.resize(nMasters,0);
					std::ifstream infile;
					const std::string fname1 = prefix + ".outdegree." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(partID);
					infile.open(fname1.c_str(), std::ifstream::in | std::ifstream::binary);	
					infile.read(reinterpret_cast<char*>(&outdegrees[0]), nMasters*sizeof(int));
					infile.close();
					std::cout<<partID<<"|finish loading"<<fname1<<std::endl<<std::flush;
				} while(0);	

				//2
				do {
					const std::string fname2 = prefix + ".scatterMap." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(partID);
					std::ifstream infile;
					infile.open(fname2.c_str(), std::ifstream::in | std::ifstream::binary);	
					std::vector<Bitmap_t>& scatterMaps = distributedGraph.scatterMaps;
					scatterMaps.resize(nParts);
					for(int i=0; i<nParts; i++){
						if(i != partID){
							scatterMaps[i].in(infile);
						}
					}
					infile.close();
					std::cout<<partID<<"|finish loading"<<fname2<<std::endl<<std::flush;
				} while(0);
				//3
				do {
					const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(partID);
					std::ifstream infile;
					infile.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
					infile.seekg(0, std::ifstream::end);
					const size_t nE = (size_t)infile.tellg()/sizeof(edge_t);
					edges.resize(nE);
					infile.seekg(0, std::ifstream::beg);
					infile.read(reinterpret_cast<char*>(&edges[0]), nE*sizeof(edge_t));
					infile.close();
					std::cout<<partID<<"|finish loading"<<fname<<"..."<<nE<<"edges"<<std::endl<<std::flush;
					//std::cout<<partID<<"sample|"<<edges[0].first<<"-->"<<edges[0].second<<std::endl<<std::flush;
				} while(0);
			}
			//
			//Internal Classes
			//
			//
			class GreedyEdgeDecision{
			private:
				std::vector<FixedBitmap_t>& srcMap;
				std::vector<FixedBitmap_t>& dstMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				const int partID;
				//Local Data Structures
				static const double epsilon = 1.0;
				std::vector<double> vScores;
			public:
				GreedyEdgeDecision(const int nParts_, const int partID_, std::vector<FixedBitmap_t>& srcMap_, std::vector<FixedBitmap_t>& dstMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),partID(partID_), srcMap(srcMap_),dstMap(dstMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}

				~GreedyEdgeDecision(){}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					size_t minedges = GRE::Util::getMin(edgeCounts);
					size_t maxedges = GRE::Util::getMax(edgeCounts);
					const double base = epsilon + (double)(maxedges - minedges);
					for(int i=0; i<nParts; i++){
						const double bal = 1+(maxedges - edgeCounts[i])/base;
						//vScores[i] = bal * ((srcMap[src].test(i)? 1.0:0.0) + (src%nParts==i? 0.5:0.0)
						//+ (dstMap[dest].test(i)? 1.0:0.0) + (dest%nParts==i? 0.5:0.0));
						vScores[i] = bal * ((srcMap[src].test(i)? 1.0:(src%nParts==i? 0.5:0.0)
						+ (dstMap[dest].test(i)? 1.0:(dest%nParts==i? 0.5:0.0))));
						//vScores[i] = bal + (srcMap[src].test(i)? 1.0:(src%nParts==i? 0.5:0.0))
						//+ (dstMap[dest].test(i)? 1.0:(dest%nParts==i? 0.5:0.0));
						//vScores[i] = bal + (srcMap[src].test(i)? 1.0:0.0)+(dstMap[dest].test(i)? 1.0:0.0);

					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					int pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					srcMap[src].set_serial(pid);
					dstMap[dest].set_serial(pid);
					return pid;
				}
			};
			class LDGEdgeDecision{
			private:
				std::vector<FixedBitmap_t>& srcMap;
				std::vector<FixedBitmap_t>& dstMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				const int partID;
				//Local Data Structures
				const long capacity;
				std::vector<double> vScores;
			public:
				LDGEdgeDecision(const int nParts_, const int partID_, const long capacity_, std::vector<FixedBitmap_t>& srcMap_, std::vector<FixedBitmap_t>& dstMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),partID(partID_), capacity(capacity_), srcMap(srcMap_),dstMap(dstMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}
				~LDGEdgeDecision(){}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					for(int i=0; i<nParts; i++){
						const double bal = fabs(1.0 - (double)edgeCounts[i]/(double)capacity);
						//vScores[i] = bal * ((srcMap[src].test(i)? 1.0:0.0) + (src%nParts==i? 0.5:0.0)
						//+ (dstMap[dest].test(i)? 1.0:0.0) + (dest%nParts==i? 0.5:0.0));
						vScores[i] = bal * ((srcMap[src].test(i)? 1.0:(src%nParts==i? 0.5:0.0)
						+ (dstMap[dest].test(i)? 1.0:(dest%nParts==i? 0.5:0.0))));
						//vScores[i] = bal + (srcMap[src].test(i)? 1.0:(src%nParts==i? 0.5:0.0))
						//+ (dstMap[dest].test(i)? 1.0:(dest%nParts==i? 0.5:0.0));
						//vScores[i] = bal + (srcMap[src].test(i)? 1.0:0.0)+(dstMap[dest].test(i)? 1.0:0.0);

					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					int pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					srcMap[src].set_serial(pid);
					dstMap[dest].set_serial(pid);
					return pid;
				}
			};


			class EdgeDecision{
			private:
				std::vector<FixedBitmap_t>& srcMap;
				std::vector<FixedBitmap_t>& dstMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				const int partID;
				//Local Data Structures
				static const double epsilon = 1.0;
				std::vector<double> vScores;
			public:
				EdgeDecision(const int nParts_, const int partID_, std::vector<FixedBitmap_t>& srcMap_, std::vector<FixedBitmap_t>& dstMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),partID(partID_),srcMap(srcMap_),dstMap(dstMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}

				~EdgeDecision(){}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					size_t minedges = GRE::Util::getMin(edgeCounts);
					size_t maxedges = GRE::Util::getMax(edgeCounts);
					const double base = epsilon + (double)(maxedges - minedges);
					for(int i=0; i<nParts; i++){
						const double bal = (maxedges - edgeCounts[i])/base;
						vScores[i] = bal + (srcMap[src].test(i)? 1.0:0.0) + (src%nParts==i? 0.5:0.0)
						+ (dstMap[dest].test(i)? 1.0:0.0) + (dest%nParts==i? 0.5:0.0);
					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					int pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					srcMap[src].set_serial(pid);
					dstMap[dest].set_serial(pid);
					return pid;
				}
			};

			class EdgeDecision1{
			private:
				std::vector<FixedBitmap_t>& vtxMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				const int partID;
				//Local Data Structures
				static const double epsilon = 1.0;
				std::vector<double> vScores;
			public:
				EdgeDecision1(const int nParts_, const int partID_, std::vector<FixedBitmap_t>& vtxMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),partID(partID_),vtxMap(vtxMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}

				~EdgeDecision1(){}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					size_t minedges = GRE::Util::getMin(edgeCounts);
					size_t maxedges = GRE::Util::getMax(edgeCounts);
					const double base = epsilon + (double)(maxedges - minedges);
					for(int i=0; i<nParts; i++){
						const double bal = (maxedges - edgeCounts[i])/base;
						vScores[i] = bal + (vtxMap[src].test(i)? 1.0:0.0) + (src%nParts==i? 0.5:0.0)
							+ (vtxMap[dest].test(i)? 1.0:0.0) + (dest%nParts==i? 0.5:0.0);
				
					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					int pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					vtxMap[src].set_serial(pid);
					vtxMap[dest].set_serial(pid);
					return pid;
				}
			};

			inline int getPartID(const vid_t vid, const int nParts){
				return vid % nParts;
			}
		};

		//
		//Partition edge lists to n parts
		//Single Process version, optimal
		class EdgePartition {
		private:
			typedef GRE::Util::Bitmap Bitmap_t;
			typedef GRE::Util::FixedBitmap<GRE::MAX_NUM_MACHINES> FixedBitmap_t;
		private:
			size_t nV;
			size_t nE;
		public:
			EdgePartition():nV(0), nE(0){}
			~EdgePartition(){}
		public:
			void divideEdges(const std::string& prefix, const int n){
				loadMetaFromFile(prefix);
				redistributeEdgeFiles(prefix, n);
			}
		private:
			void loadMetaFromFile(const std::string& prefix){
				std::string filename = prefix + ".meta";
				std::ifstream infile(filename.c_str());
				infile >> nV;
				infile >> nE;
			}
			//three output files:
			//1: edgelist.num.id
			//2: outdegree.num.id
			//3: scattermap.num.id
			void redistributeEdgeFiles(const std::string& prefix, const int nParts){
				//open stream
				//read-in
				const size_t nVertices = nV;
				const size_t nEdges = nE;
			
				static size_t bufSize = 0x800000L;
				//stream buf 3: edgelist files
				std::vector<edge_t> edgeLists_buf[nParts];
				for(int i = 0; i < nParts; i++){
					edgeLists_buf[i].reserve(bufSize);
				}

				std::vector<FixedBitmap_t> srcMap;
				std::vector<FixedBitmap_t> dstMap;
				std::vector<size_t> edgeCounts;

				srcMap.resize(nVertices);
				dstMap.resize(nVertices);
				edgeCounts.resize(nParts, 0);
			
				EdgeDecision edgeDecision(nParts, srcMap, dstMap, edgeCounts);

				std::vector<int> outdegree_buf;
				std::vector<Bitmap_t> scatterMap_buf;
				scatterMap_buf.resize(nParts);
				for(int p = 0; p < nParts; p++){
					//
					const size_t nMasters = nVertices/nParts + (nVertices%nParts > p? 1:0);
					outdegree_buf.resize(nMasters, 0);
					for(int i = 0; i < nParts; i++){
						scatterMap_buf[i].resize(nMasters);
					}
					//
					for(int i=p; ;i+=nParts)
					{
						std::string filename(prefix);
						filename = filename + ".edges."+ GRE::Util::toString(i);
						std::ifstream infile(filename.c_str());
						if(!infile) {
							std::cerr<<p<<" : Ok! Finish loading edgelist files."<<std::endl;
							break;
						}
						std::cout<<"Loading edges from"<<filename<<"..."<<std::endl<<std::flush;
						//SNAP and TSV
						while(infile.good() && !infile.eof()){
							std::string line;
							std::getline(infile, line);
							if(line.empty() || line[0]=='#') continue;
							if(infile.fail()) break;
							char* tmpptr;
							const vid_t src = strtoul(line.c_str(), &tmpptr, 10);
							const vid_t dest = strtoul(tmpptr, NULL, 10);
							int epid = edgeDecision.getPartID(src, dest);			
							const edge_t edge(src, dest);
							edgeLists_buf[epid].push_back(edge);
							if(edgeLists_buf[epid].size()==bufSize){
								const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(epid);
								std::ofstream outfile;
								outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
								outfile.write(reinterpret_cast<char*>(&edgeLists_buf[epid][0]), bufSize*sizeof(edge_t));
								outfile.close();
								edgeLists_buf[epid].clear();
							}
							//
							lvid_t lsrc = src/nParts;
							outdegree_buf[lsrc]++;
							scatterMap_buf[epid].set_serial(lsrc);
						}//end-while-reading
					}
					do {
						std::ofstream outfile;
						//1
						const std::string fname1 = prefix + ".outdegree." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
						outfile.open(fname1.c_str(), std::ofstream::out | std::ofstream::binary);	
						outfile.write(reinterpret_cast<char*>(&outdegree_buf[0]), nMasters*sizeof(int));
						outfile.close();
						outdegree_buf.clear();
						//2
						const std::string fname2 = prefix + ".scatterMap." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
						outfile.open(fname2.c_str(), std::ofstream::out | std::ofstream::binary);	
						for(int j=0; j<nParts; j++){
							if(j != p)
								scatterMap_buf[j].out(outfile);
							scatterMap_buf[j].clear();
						}
						outfile.close();
					} while(0);
				}//end-for-reading-all
				//flush
				for(int p=0; p<nParts; p++){
					if(edgeLists_buf[p].size()==0) continue;
					const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
					std::ofstream outfile;
					outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
					outfile.write(reinterpret_cast<char*>(&edgeLists_buf[p][0]), edgeLists_buf[p].size()*sizeof(edge_t));
					outfile.close();
					edgeLists_buf[p].clear();
				}		
				//optional sync
				srcMap.clear();
				dstMap.clear();
				edgeCounts.clear();
			}
			//
			//Internal Class
			//


			class EdgeDecision{
			private:
				std::vector<FixedBitmap_t>& srcMap;
				std::vector<FixedBitmap_t>& dstMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				//Local Data Structures
				static const double epsilon = 1.0;
				std::vector<double> vScores;
			public:
				EdgeDecision(const int nParts_, std::vector<FixedBitmap_t>& srcMap_, std::vector<FixedBitmap_t>& dstMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),srcMap(srcMap_),dstMap(dstMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}
				~EdgeDecision(){
					vScores.clear();
				}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					size_t minedges = GRE::Util::getMin(edgeCounts);
					size_t maxedges = GRE::Util::getMax(edgeCounts);
					const double base = epsilon + (double)(maxedges - minedges);
					for(int i=0; i<nParts; i++){
						const double bal = (maxedges - edgeCounts[i])/base;
						vScores[i] = bal + (srcMap[src].test(i)? 1.0:0.0) + (dstMap[dest].test(i)? 1.0:0.0);
					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					const int pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					srcMap[src].set_serial(pid);
					dstMap[dest].set_serial(pid);
					return pid;
				}
			};
		};//end edgePar

		//GRE-P
		//Partition edge lists to n parts
		//Multiple Process version, not optimal
		class parEdgePartition{
		public:
			typedef GRE::COMM::Underlying Comm_t;
		private:
			typedef GRE::Util::Bitmap Bitmap_t;
			typedef GRE::Util::FixedBitmap<GRE::MAX_NUM_MACHINES> FixedBitmap_t;
		private:
			Comm_t& comm;
			size_t nV;
			size_t nE;
		public:
			parEdgePartition(Comm_t& _comm):comm(_comm),nV(0),nE(0){}
			~parEdgePartition(){}
		public:
			void divideEdges(const std::string& prefix){
				loadMetaFromFile(prefix);
				redistributeEdgeFiles(prefix);
			}
		private:
			void loadMetaFromFile(const std::string& prefix){
				std::string filename = prefix + ".meta";
				std::ifstream infile(filename.c_str());
				infile >> nV;
				infile >> nE;
			}
			//three output files:
			//1: edgelist.num.id
			//2: outdegree.num.id
			//3: scattermap.num.id
			void redistributeEdgeFiles(const std::string& prefix){
				//open stream
				//read-in
				const size_t nVertices = nV;
				const size_t nEdges = nE;
				const int nParts = comm.size();
				const int partID = comm.rank();
				static size_t bufSize = 0x800000L;

				//stream buf 3: edgelist files
				std::vector<edge_t> edgeLists_buf[nParts];
				for(int i = 0; i < nParts; i++){
					edgeLists_buf[i].reserve(bufSize);
				}

				std::vector<FixedBitmap_t> srcMap;
				std::vector<FixedBitmap_t> dstMap;
				std::vector<size_t> edgeCounts;

				srcMap.resize(nVertices);
				dstMap.resize(nVertices);
				edgeCounts.resize(nParts, 0);
			
				EdgeDecision edgeDecision(nParts, partID, srcMap, dstMap, edgeCounts);

				std::vector<int> outdegree_buf;
				std::vector<Bitmap_t> scatterMap_buf;
				scatterMap_buf.resize(nParts);
		
				const int p = comm.rank();
				do{
					const size_t nMasters = nVertices/nParts + (nVertices%nParts > p? 1:0);
					outdegree_buf.resize(nMasters, 0);
					for(int i = 0; i < nParts; i++){
						scatterMap_buf[i].resize(nMasters);
					}
					//
					for(int i=p; ;i+=nParts)
					{
						std::string filename(prefix);
						filename = filename + ".edges."+ GRE::Util::toString(i);
						std::ifstream infile(filename.c_str());
						if(!infile) {
							std::cerr<<p<<" : Ok! Finish loading edgelist files."<<std::endl;
							break;
						}
						std::cout<<"Loading edges from"<<filename<<"..."<<std::endl<<std::flush;
						//SNAP and TSV
						while(infile.good() && !infile.eof()){
							std::string line;
							std::getline(infile, line);
							if(line.empty() || line[0]=='#') continue;
							if(infile.fail()) break;
							char* tmpptr;
							const vid_t src = strtoul(line.c_str(), &tmpptr, 10);
							const vid_t dest = strtoul(tmpptr, NULL, 10);
							int epid = edgeDecision.getPartID(src, dest);			
							const edge_t edge(src, dest);
							edgeLists_buf[epid].push_back(edge);
							if(edgeLists_buf[epid].size()==bufSize){
								const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(epid) + "." + GRE::Util::toString(p);//Potential a bug when no such a file exists!
								std::ofstream outfile;
								outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
								outfile.write(reinterpret_cast<char*>(&edgeLists_buf[epid][0]), bufSize*sizeof(edge_t));
								outfile.close();
								edgeLists_buf[epid].clear();
							}
							//
							lvid_t lsrc = src/nParts;
							outdegree_buf[lsrc]++;
							scatterMap_buf[epid].set_serial(lsrc);
						}//end-while-reading
					}
					do {
						std::ofstream outfile;
						//1
						const std::string fname1 = prefix + ".outdegree." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
						outfile.open(fname1.c_str(), std::ofstream::out | std::ofstream::binary);	
						outfile.write(reinterpret_cast<char*>(&outdegree_buf[0]), nMasters*sizeof(int));
						outfile.close();
						outdegree_buf.clear();
						//2
						const std::string fname2 = prefix + ".scatterMap." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
						outfile.open(fname2.c_str(), std::ofstream::out | std::ofstream::binary);	
						for(int j=0; j<nParts; j++){
							if(j != p)
								scatterMap_buf[j].out(outfile);
							scatterMap_buf[j].clear();
						}
						outfile.close();
					} while(0);
				} while(0);//end-for-reading-all
				//flush
				for(int i=0; i<nParts; i++){
					if(edgeLists_buf[i].size()==0) continue;
					const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(i) 
						+ "." + GRE::Util::toString(p);//Potential a bug when no such a file exists!
					std::ofstream outfile;
					outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
					outfile.write(reinterpret_cast<char*>(&edgeLists_buf[i][0]), edgeLists_buf[i].size()*sizeof(edge_t));
					outfile.flush();
					outfile.close();
					edgeLists_buf[i].clear();
				}		
				//sync
				comm.barrier();
				//merge
				do{
					std::vector<edge_t> edges;
					const size_t bufSizeX = 8*bufSize;
					edges.resize(bufSizeX);
					const std::string fname = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p);
					std::ofstream outfile;
					outfile.open(fname.c_str(), std::ofstream::out | std::ofstream::app | std::ofstream::binary);
					for(int i=0; i<nParts; i++){
						const std::string fname1 = prefix + ".part." + GRE::Util::toString(nParts) + "." + GRE::Util::toString(p) 
							+ "." + GRE::Util::toString(i);//Potential a bug when no such a file exists!
						std::ifstream infile;
						infile.open(fname1.c_str(), std::ifstream::in|std::ifstream::binary);
						infile.seekg(0, std::ifstream::end);
						const size_t size = infile.tellg();
						const int n = size / (bufSizeX*sizeof(edge_t));
						const size_t left = size % (bufSizeX*sizeof(edge_t));
						infile.seekg(0, std::ifstream::beg);
						for(int j = 0; j < n; j++){
							infile.read(reinterpret_cast<char*>(&edges[0]), bufSizeX*sizeof(edge_t));
							outfile.write(reinterpret_cast<char*>(&edges[0]), bufSizeX*sizeof(edge_t));
						}
						if(left > 0){
							infile.read(reinterpret_cast<char*>(&edges[0]), left*sizeof(edge_t));
							outfile.write(reinterpret_cast<char*>(&edges[0]), left*sizeof(edge_t));
						}
						infile.close();
					}		
					outfile.close();
				} while(0);
				//optional sync
				srcMap.clear();
				dstMap.clear();
				edgeCounts.clear();
			}
			//
			//Internal Class
			//
			class EdgeDecision{
			private:
				std::vector<FixedBitmap_t>& srcMap;
				std::vector<FixedBitmap_t>& dstMap;
				std::vector<size_t>& edgeCounts;
				const int nParts;
				const int partID;
				//Local Data Structures
				static const double epsilon = 1.0;
				std::vector<double> vScores;
			public:
				EdgeDecision(const int nParts_, const int partID_, std::vector<FixedBitmap_t>& srcMap_, std::vector<FixedBitmap_t>& dstMap_, std::vector<size_t>& edgeCounts_)
					:nParts(nParts_),partID(partID_),srcMap(srcMap_),dstMap(dstMap_),edgeCounts(edgeCounts_)
				{
					vScores.resize(nParts, 0.0);
				}
				~EdgeDecision(){
					vScores.clear();
				}
			public:
				inline int getPartID(const vid_t src, const vid_t dest)
				{
					int pid=getPartID(src, nParts);
					if(pid==getPartID(dest, nParts)) {
						edgeCounts[pid]++;
						//srcMap[src].set_serial(pid);
						//dstMap[dest].set_serial(pid);
						return pid;
					}
					size_t minedges = GRE::Util::getMin(edgeCounts);
					size_t maxedges = GRE::Util::getMax(edgeCounts);
					const double base = epsilon + (double)(maxedges - minedges);
					for(int i=0; i<nParts; i++){
						const double bal = (maxedges - edgeCounts[i])/base;
						vScores[i] = bal + (srcMap[src].test(i)? 1.0:0.0) + (dstMap[dest].test(i)? 1.0:0.0);
					}
					const double maxScore = GRE::Util::getMax(vScores);
					int stack[nParts];
					int stackTop = 0;
					for(int i=0; i<nParts; i++){
						if(fabs(maxScore - vScores[i]) < 1e-4){
							stack[stackTop++]=i;
						}
					}
					pid = stack[(src+dest)%stackTop];//make the choice random.
					edgeCounts[pid]++;
					srcMap[src].set_serial(pid);
					dstMap[dest].set_serial(pid);
					return pid;
				}
			};
		};//end edgePar
	}//end Namespace-Graph
}
#endif
