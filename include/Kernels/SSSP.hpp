#ifndef _GRE_APP_SSSP_HPP
#define _GRE_APP_SSSP_HPP

#include "Types.hpp"
#include "Util/vLock.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"

typedef uint32_t dist_t;
typedef dist_t vertex_data_t;
typedef dist_t scatter_data_t;
typedef dist_t combine_data_t;

typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> SSSP_Base;

class SSSP: public SSSP_Base{
public:
	SSSP(DistributedGraph_t& dg):SSSP_Base(dg), states(NULL),lastStates(NULL), edgeWgt(NULL){}
	~SSSP(){}

	void alloc(){
		allocScatterData(lastStates);
		allocCombineData(states, std::numeric_limits<dist_t>::max());
		allocVertexData(states);
	}
	void free(){
		freeScatterData(lastStates);
		freeCombineData(states);
		freeVertexData(states);
	}
	void loadVertexProperty(const std::string prefix){
		//
	}
	void unloadVertexProperty(){
		//
	}
	void loadEdgeProperty(const std::string filename){
		graph.loadEdgeProperty(edgeWgt, filename);
	}
	void loadEdgeProperty(){
		graph.loadEdgeProperty(edgeWgt);
	}
	void unloadEdgeProperty(){
		graph.unloadEdgeProperty(edgeWgt);
	}
private:
	vertex_data_t* states;
	scatter_data_t* lastStates;
	dist_t* edgeWgt;
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		states[lv] = 0;
		lastStates[lv] = states[lv];//same with states.
		return true;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		const leid_t leid = graph.refLocalGraph().leid(it);
		const dist_t dist = lastStates[lsrc] + edgeWgt[leid];
		scatter_data_t msg_data(dist);
		ctx.engine->sendMessage(ctx, ldst, msg_data);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		if(msg_data < states[ldst]){
			lock.spinLock();
			if(msg_data < states[ldst]){
				states[ldst] = msg_data;
				lock.unLock();
				return true;
			}
			lock.unLock();
		}
		return false;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		lastStates[lv] = states[lv];
		return true;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
	}
	inline void reset_scatter(const lvid_t lsrc){
	}
}; 
#endif
