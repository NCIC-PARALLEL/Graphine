#ifndef _GRE_APP_BC_HPP
#define _GRE_APP_BC_HPP

#include "Types.hpp"
#include "Util/vLock.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"

//vertex_data_t for BFS
struct vtxData{
	int level;
	int numPaths;	
	vtxData():level(-1), numPaths(0){}
	~vtxData(){}
};
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<int, int, vtxData, DistributedGraph_t> Forward_Base;
typedef	GRE::VertexCode<double, double, double, DistributedGraph_t> Backward_Base;

class Forward: public Forward_Base{
public:
	Forward(DistributedGraph_t& dg, vtxData*& _state, int _currLevel=1):Forward_Base(dg), state(_state), numPaths(NULL), sum(NULL),
		currLevel(_currLevel){}

	~Forward(){}

	void alloc(){
		allocScatterData(numPaths, 0);
		allocCombineData(sum, 0);
		allocVertexData(state);
	}
	void free(){
		freeScatterData(numPaths);
		freeCombineData(sum);
		//freeVertexData(state);
	}
	template<typename dc_t>
	void reduce(dc_t& dc){
		currLevel++;
	}
	void setCurrLevel(int l){currLevel=l;}
	int getCurrLevel(){return currLevel;}
private:
	//vertex_data
	vtxData*& state;
	//scatter data
	int* numPaths;
	//combine data
	int* sum;
	//other
	int currLevel;
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		//update my state
		state[lv].numPaths=1;
		state[lv].level=0;
		//update scatter data
		numPaths[lv]=state[lv].numPaths;
		return true;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, numPaths[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		sum[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(state[lv].level<0){//unvisited
			//update my state
			state[lv].level=currLevel;
			state[lv].numPaths=sum[lv];
			//update scatter data
			numPaths[lv]=state[lv].numPaths;
			return true;
		}
		return false;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){}
}; 

class Backward: public Backward_Base{
public:
	Backward(DistributedGraph_t& dg, vtxData* _bfsState, double*& _dependencies, int _currLevel=0):
		Backward_Base(dg), bfsState(_bfsState), dependencies(_dependencies), 
		pathwgt(NULL), sum(NULL), currLevel(_currLevel){}
	~Backward(){}

	void setCurrLevel(const int l){ currLevel=l;}
	int getCurrLevel(){ return currLevel;}
	void aggregate(){}

	void alloc(){
		assert(bfsState!=NULL);
		allocScatterData(pathwgt, 0);
		allocCombineData(sum, 0);
		allocVertexData(dependencies, 0.0);
	}
	void free(){
		freeScatterData(pathwgt);
		freeCombineData(sum);
		//freeVertexData(dependencies);
	}
private:
	//external
	vtxData* bfsState;
	//vertex_data
	double*& dependencies;
	//scatter
	double* pathwgt;
	//combine
	double* sum;
	//other
	int currLevel;
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		if(bfsState[lv].level==currLevel){
			pathwgt[lv]=(1.0+dependencies[lv])/bfsState[lv].numPaths;
			return true;
		} else
			return false;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, pathwgt[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		sum[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(bfsState[lv].level==currLevel-1){
			dependencies[lv]=sum[lv]*bfsState[lv].numPaths;
			//pathwgt[lv]=(1.0+dependencies[lv])/numPaths[lv];
			//return true;
		}
		return false;
	}
	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
		sum[ldst]=0.0;
	}
}; 

#endif
