#ifndef _GRE_APP_CC_HPP
#define _GRE_APP_CC_HPP

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
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> ConnectedComponents_Base;

class ConnectedComponents: public ConnectedComponents_Base{
public:
	ConnectedComponents(DistributedGraph_t& dg): ConnectedComponents_Base(dg), cc_scatter(NULL), cc_combine(NULL), cc(NULL){}
	~ConnectedComponents(){}
	void alloc(){
		allocScatterData(cc_scatter);
		allocCombineData(cc_combine, std::numeric_limits<vid_t>::max());
		allocVertexData(cc, std::numeric_limits<vid_t>::max());
	}
	void free(){
		freeScatterData(cc_scatter);
		freeCombineData(cc_combine);
		freeVertexData(cc);
	}
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		cc[lv] = graph.get_vid(lv);
		cc_scatter[lv] = cc[lv];
		cc_combine[lv] = cc[lv];//optioanl optimization: speed up.
		return true;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, cc_scatter[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		if(cc_combine[ldst]>msg_data){
			lock.spinLock();
			if(cc_combine[ldst]>msg_data){
				cc_combine[ldst]=msg_data;
				lock.unLock();
				return true;
			}
			lock.unLock();
		}
		return false;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(cc_combine[lv]>=cc[lv]) return false;
		cc[lv] = cc_combine[lv];
		cc_scatter[lv] = cc[lv];
		return true;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){}

	/////////////Extra///////////////
	long get_num_components(){
		const size_t n = graph.get_num_masters();
		long nComponents = 0;
		/*
 		#pragma omp parallel for reduction(+:nComponents)
		for(lv = 0; lv < n; lv++){
			if(cc[lv] == graph.get_vid(lv))
				nComponents++;
		}
		*/
		#pragma omp parallel
		{
			long myN = 0;
			lvid_t lv;
 			#pragma omp for
			for(lv = 0; lv < n; lv++){
				if(cc[lv] == graph.get_vid(lv))
					myN++;
			}
			AtomicFetchAdd(&nComponents, myN);
		}
		return nComponents; 
	}

	template <typename DC_T>
	void reduce(DC_T& dc){
		long nComponents = get_num_components();
		long nC = graph.get_num_combiners();
		long nS = graph.get_num_scatters();
		long nE = graph.get_num_edges();
		long allNComponents = 0;
		long allNC = 0;
		long allNS = 0;
		long allNE = 0;
		const int root = 0;
		dc.refUnderlying().reduce(&nComponents, &allNComponents, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nC, &allNC, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nS, &allNS, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nE, &allNE, 1, MPI_LONG, MPI_SUM, root);
		if(dc.procID() == root){
			std::cout << "Number of components: "<< allNComponents << std::endl;
			std::cout << "Number of combiners: "<< allNC << std::endl;
			std::cout << "Number of scatters: "<< allNS << std::endl;
			std::cout << "Number of agents: "<< allNC+allNS << std::endl;
			std::cout << "Number of edges: "<< allNE << std::endl;
			std::cout << "Value of equevalent cut: "<< (double)(allNC+allNS)/(double)(allNE) << std::endl;
		}
	}
private:
	vertex_data_t* cc;
	scatter_data_t* cc_scatter;
	combine_data_t* cc_combine;
};

#endif
