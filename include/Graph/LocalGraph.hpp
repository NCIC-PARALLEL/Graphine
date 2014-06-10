#ifndef GRE_Graph_LocalGraph_HPP
#define GRE_Graph_LocalGraph_HPP

#include <iostream>

#include "Types.hpp"
#include "Graph/GraphIngress.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Util/Atomic.hpp"
//#include "Util/IOArchive.hpp"
#include <vector>
#include <omp.h>

namespace GRE{
	namespace Graph{
		//namespace util = GRE::Util;
		//
		//Graph is used as an Index.
		//
		template<typename T> class DistributedGraph;
		template<typename T> class GraphIngress;

		class GraphInCSR {
		public:
			friend class DistributedGraph<GraphInCSR>;
			friend class GraphIngress<GraphInCSR>;
			typedef std::vector<lidx_t> EdgeIndex_t;
			typedef std::vector<lvid_t> EdgeList_t;
		private:
			//ready flag
			bool isReady;
			bool isrReady;
			//meta-data
			size_t nVertices;
			size_t nEdges;
			//graph structure:Use CSR format
			EdgeIndex_t edgeIndex;
			EdgeList_t edgeList;
			//reserve CSR
			EdgeIndex_t redgeIndex;
			EdgeList_t redgeList;
			EdgeIndex_t redgeID;//internal, never use any other where.

		public:
			GraphInCSR():isReady(false), isrReady(false), nVertices(0), nEdges(0), edgeIndex(), edgeList(), redgeIndex(), redgeList(){}
			~GraphInCSR(){}
		public:
			/*
			void load(util::iarchive& iarc)
			{
				iarc >> nVertices
				util::iread(iarc, edgeIndex);
				iarc >> nEdges;
				util::iread(iarc, edgeList);
				isReady = true;
			}
			void save(util::oarchive& oarc)
			{
				oarc << nVertices
				util::iwrite(edgeIndex);
				oarc << nEdges;
				util::iwrite(edgeList);
			}
			*/
			void build(const std::vector<ledge_t>& edges){}
			void build(const std::vector<ledge_t>& edges, const size_t nV){
				const size_t nE = edges.size();
				edgeIndex.resize(nV+1, 0);
				edgeList.resize(nE);

				#pragma omp parallel for
				for(lidx_t i=0; i<nE; i++){
					const lvid_t lsrc = edges[i].first;
					AtomicFetchAdd(&edgeIndex[lsrc], 1);
				}
				lidx_t off = edgeIndex[0];
				edgeIndex[0] = 0;
				for(lvid_t i=1; i<nV+1; i++){
					const lidx_t off1 = off;
					off = edgeIndex[i];
					edgeIndex[i] = off1+edgeIndex[i-1];
				}

				#pragma omp parallel 
				{	
					const int nThreads = omp_get_num_threads();
					const int tid = omp_get_thread_num();
					for(lidx_t i=0; i<nE; i++){
						const lvid_t lsrc = edges[i].first;
						if(lsrc % nThreads == tid){
							edgeList[edgeIndex[lsrc]++]=edges[i].second;
						}
					}
				}
				for(lvid_t i=nV; i>0; i--){
					edgeIndex[i] = edgeIndex[i-1];
				}
				edgeIndex[0] = 0;
				if(nVertices == 0) nVertices = nV;
				assert(nVertices == nV);
				if(nEdges == 0) nEdges = nE;
				assert(nEdges == nE);
				//
				//std::cout<<"# of Local Vertices is "<<nV<<"."<<std::endl;
				//std::cout<<"# of Local Edges is "<<nE<<"."<<std::endl;
				isReady = true;
			}
			void rbuild(const std::vector<ledge_t>& edges){
				//TODO
			}

			void rbuild(const std::vector<ledge_t>& edges, const size_t nV){
				const size_t nE = edges.size();
				redgeIndex.resize(nV+1, 0);
				redgeList.resize(nE);
				redgeID.resize(nE);

				#pragma omp parallel for
				for(lidx_t i=0; i<nE; i++){
					const lvid_t lsrc = edges[i].second;
					AtomicFetchAdd(&redgeIndex[lsrc], 1);
				}
				lidx_t off = redgeIndex[0];
				redgeIndex[0] = 0;
				for(lvid_t i=1; i<nV+1; i++){
					const lidx_t off1 = off;
					off = redgeIndex[i];
					redgeIndex[i] = off1+redgeIndex[i-1];
				}

				#pragma omp parallel 
				{	
					const int nThreads = omp_get_num_threads();
					const int tid = omp_get_thread_num();
					for(lidx_t i=0; i<nE; i++){
						const lvid_t lsrc = edges[i].second;
						if(lsrc % nThreads == tid){
							redgeList[redgeIndex[lsrc]]=edges[i].first;
							redgeID[redgeIndex[lsrc]]=i;
							redgeIndex[lsrc]++;
						}
					}
				}
				for(lvid_t i=nV; i>0; i--){
					redgeIndex[i] = redgeIndex[i-1];
				}
				redgeIndex[0] = 0;
				if(nVertices == 0) nVertices = nV;
				assert(nVertices == nV);
				if(nEdges == 0) nEdges = nE;
				assert(nEdges == nE);
				//
				//std::cout<<"# of Local Vertices is "<<nV<<"."<<std::endl;
				//std::cout<<"# of Local Edges is "<<nE<<"."<<std::endl;
				isrReady = true;
			}

		public:
			//edgeIterator
			typedef EdgeList_t::iterator edge_iterator;
			inline edge_iterator edge_begin(){
				return edgeList.begin();
			}
			inline edge_iterator edge_begin(lvid_t lvid){
				return edgeList.begin()+edgeIndex[lvid];
			}
			inline edge_iterator edge_end(){
				return edgeList.end();
			}
			inline edge_iterator edge_end(lvid_t lvid){
				return edgeList.begin()+edgeIndex[lvid+1];
			}
			inline leid_t leid(const edge_iterator& it){
				const lidx_t idx = (it - this->edge_begin())/sizeof(lidx_t);
				return idx;
			}
			//
			//reverse edgeIterator
			typedef EdgeList_t::iterator redge_iterator;
			inline redge_iterator redge_begin(){
				return redgeList.begin();
			}
			inline redge_iterator redge_begin(lvid_t lvid){
				return redgeList.begin()+edgeIndex[lvid];
			}
			inline redge_iterator redge_end(){
				return redgeList.end();
			}
			inline redge_iterator redge_end(lvid_t lvid){
				return redgeList.begin()+edgeIndex[lvid+1];
			}
			inline leid_t rleid(const redge_iterator& it){
				const lidx_t idx = redgeID[(it - this->redge_begin())/sizeof(lidx_t)];
				return idx;
			}

			inline size_t get_num_edges(){
				return nEdges;
			} 
		};
		//
		//Property Graph
		//TODO:not used in GRE.
		//
		template <typename VertexProperty, typename EdgeProperty>
		class PropertyGraphInCSR{
		private:
			PropertyGraphInCSR(){/*not defined*/}
		};
	}//end NameSpce-Graph
}
#endif
