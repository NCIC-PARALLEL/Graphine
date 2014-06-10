#ifndef _GRE_APP_CountIn_HPP
#define _GRE_APP_CountIn_HPP

#include "Types.hpp"
#include "Util/vLock.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"

typedef int vertex_data_t;
typedef int scatter_data_t;
typedef int combine_data_t;

typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> CountIn_Base;

class CountIn: public CountIn_Base{
public:
	CountIn(DistributedGraph_t& dg): CountIn_Base(dg), cn(NULL), cn_scatter(NULL), cn_combine(NULL){}
	~CountIn(){}
	void alloc(){
		allocScatterData(cn_scatter);
		allocCombineData(cn_combine, 0);
		allocVertexData(cn);
	}
	void free(){
		freeScatterData(cn_scatter);
		freeCombineData(cn_combine);
		freeVertexData(cn);
	}
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		cn[lv] = 0;
		cn_scatter[lv] = 1;
		return true;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, cn_scatter[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		cn_combine[ldst]+=msg_data;
		//cn_combine[ldst]++;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		cn[lv]=cn_combine[lv];
		return false;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
		//do nothing
	}

	/////////////Extra///////////////
	long get_num_outEdges(){
		long nOut = 0;
		vid_t n = graph.get_num_masters();
		lvid_t lv;
 		#pragma omp parallel for reduction(+:nOut)
		for(lv = 0; lv < n; lv++){
			nOut+=graph.getOutdegree(lv);
		}
		return nOut; 
	}

	long get_num_inEdges(){
		long nIn = 0;
		vid_t n = graph.get_num_masters();
		lvid_t lv;
 		#pragma omp parallel for reduction(+:nIn)
		for(lv = 0; lv < n; lv++){
			nIn+=cn[lv];
		}
		return nIn; 
	}
	template <typename DC_T>
	void output_statistics(DC_T& dc){
		long nInEdges = get_num_inEdges();	
		long nOutEdges = get_num_outEdges();	
		long allIn = 0;
		long allOut = 0;
		const int root = 0;
		dc.refUnderlying().reduce(&nInEdges, &allIn, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nOutEdges, &allOut, 1, MPI_LONG, MPI_SUM, root);
		if(dc.procID() == root){
			std::cout << "Number of in edges: "<< allIn << std::endl;
			std::cout << "Number of out edges: "<< allOut << std::endl;
		}
	}
private:
	vertex_data_t* cn;
	scatter_data_t* cn_scatter;
	combine_data_t* cn_combine;
};

#endif
