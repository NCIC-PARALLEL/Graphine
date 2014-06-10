#ifndef _GRE_APP_BFS_HPP
#define _GRE_APP_BFS_HPP

#include "Types.hpp"
#include "Util/vLock.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"

typedef vid_t vertex_data_t;
typedef vid_t scatter_data_t;
typedef vid_t combine_data_t;

typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> BFS_Base;

class BFS: public BFS_Base{
public:
	BFS(DistributedGraph_t& dg): BFS_Base(dg), pred(NULL), pred_scatter(NULL), pred_combine(NULL){}
	~BFS(){}
	void alloc(){
		allocScatterData(pred_scatter);
		allocCombineData(pred_combine, -1);
		allocVertexData(pred, -1);
	}
	void free(){
		freeScatterData(pred_scatter);
		freeCombineData(pred_combine);
		freeVertexData(pred);
	}
	///////////////////////////////Init//////////////////////////////////
	//Source Vertex
	inline bool init(const lvid_t lv){
		pred[lv] = graph.get_vid(lv);//root!
		pred_scatter[lv] = pred[lv];
		//if(graph.getOutdegree(lv)>0){
			return true;
		//} else {
		//	return false;
		//}
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, pred_scatter[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		if(pred_combine[ldst]==-1){
			lock.spinLock();
			if(pred_combine[ldst]==-1){
				pred_combine[ldst]=msg_data;
				lock.unLock();
				return true;
			}
			lock.unLock();
		}
		return false;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(pred[lv]==-1){
			pred[lv] = pred_combine[lv];
			pred_scatter[lv] = pred[lv];
			return true;
		} else return false;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
		//do nothing
	}
	/////////////Extra///////////////
	template <typename DC_T>
	void output_stat(DC_T& dc){
	}
private:
	vertex_data_t* pred;
	scatter_data_t* pred_combine;
	scatter_data_t* pred_scatter;
};

#endif
