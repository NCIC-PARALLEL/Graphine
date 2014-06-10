#ifndef _GRE_APP_PAGERANK_HPP
#define _GRE_APP_PAGERANK_HPP

#include "Util/vLock.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"
#include <iomanip>

typedef double vertex_data_t;
typedef double scatter_data_t;
typedef double combine_data_t;

typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> PageRank_Base;

class PageRank: public PageRank_Base{
private:
	static const double dampingFactor = 0.15;
	static const double err = 1e-5;
public:
	PageRank(DistributedGraph_t& dg):PageRank_Base(dg), pr_scatter(NULL), pr_combine(NULL), pr(NULL){}
	~PageRank(){}

	void alloc(){
		allocScatterData(pr_scatter);
		allocCombineData(pr_combine, 0.0);
		allocVertexData(pr);
	}
	void free(){
		freeScatterData(pr_scatter);
		freeCombineData(pr_combine);
		freeVertexData(pr);
	}
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		pr[lv] = dampingFactor;
		//pr_combine[lv] = 0.0;
		if(graph.getOutdegree(lv) > 0){
			pr_scatter[lv] = pr[lv]/(double)graph.getOutdegree(lv);
			return true;
		} else {
			return false;
		}
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, pr_scatter[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		pr_combine[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data){
		pr_combine[ldst]+=msg_data;
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		pr[lv] = dampingFactor + (1.0-dampingFactor)*pr_combine[lv];
		pr_combine[lv] = 0.0;
		return false;
	}

	///////////////////////////////Assert-to-halt//////////////////////////////////
	inline bool assert_to_halt(const lvid_t lv){
		return false;
	}

private:
	vertex_data_t* pr;
	combine_data_t* pr_combine;
	scatter_data_t* pr_scatter;

	///////////////////////////////Misc//////////////////////////////////
public:
	inline void reset_combiner(const lvid_t ldst){
		pr_combine[ldst] = 0.0;
	}
	double get_pr_sum(){
		const size_t n = graph.get_num_masters();
		double sum = 0.0;
		lvid_t lv;
		#pragma omp parallel for reduction(+:sum)
		for(lv = 0; lv < n; lv++){
			sum+=pr[lv];
		}
		return sum;
	}
	template <typename DC_T>
	void output_statistics(DC_T& dc){
		double isum = get_pr_sum();
		long nC = graph.get_num_combiners();
		long nS = graph.get_num_scatters();
		long nE = graph.get_num_edges();
		long allNC = 0;
		long allNS = 0;
		long allNE = 0;
		double sum = 0;
		const int root = 0;
		dc.refUnderlying().reduce(&isum, &sum, 1, MPI_DOUBLE, MPI_SUM, root);
		dc.refUnderlying().reduce(&nC, &allNC, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nS, &allNS, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nE, &allNE, 1, MPI_LONG, MPI_SUM, root);
		if(dc.procID() == root){
			std::cout << std::setprecision(15) << "Sum of pr: "<< sum << std::endl << std::flush;
			std::cout << "Number of combiners: "<< allNC << std::endl << std::flush;
			std::cout << "Number of scatters: "<< allNS << std::endl << std::flush;
			std::cout << "Number of agents: "<< allNC+allNS << std::endl << std::flush;
			std::cout << "Number of edges: "<< allNE << std::endl << std::flush;
			std::cout << "Value of equevalent cut: "<< (double)(allNC+allNS)/(double)(allNE) << std::endl << std::flush;
		}
	}

}; 


class PageRank_delta: public PageRank_Base{
private:
	static const double dampingFactor = 0.15;
	static const double err = 1e-5;
public:
	PageRank_delta(DistributedGraph_t& dg):PageRank_Base(dg), pr_scatter(NULL), pr_combine(NULL), pr(NULL){}
	~PageRank_delta(){}

	void alloc(){
		allocScatterData(pr_scatter);
		allocCombineData(pr_combine, 0.0);
		allocVertexData(pr);
	}
	void free(){
		freeScatterData(pr_scatter);
		freeCombineData(pr_combine);
		freeVertexData(pr);
	}
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		pr[lv] = dampingFactor;
		//pr_combine[lv] = 0.0;
		if(graph.getOutdegree(lv)>0){
			pr_scatter[lv] = pr[lv]/(double)graph.getOutdegree(lv);
			return true;
		} else {
			return false;
		}
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, pr_scatter[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		pr_combine[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		//compute delta
		const double delta = (1.0-dampingFactor)*pr_combine[lv];
		//update vertex state
		pr[lv] += delta;
		pr_combine[lv] = 0.0;
		//filter and activate if true
		if(delta > err && graph.getOutdegree(lv)>0){
		//if(graph.getOutdegree(lv)>0){
			pr_scatter[lv] = delta/(double)graph.getOutdegree(lv); 
			return true;
		}
		return false;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
		pr_combine[ldst] = 0.0;
	}
	double get_pr_sum(){
		const size_t n = graph.get_num_masters();
		double sum = 0.0;
		lvid_t lv;
		#pragma omp parallel for reduction(+:sum)
		for(lv = 0; lv < n; lv++){
			sum+=pr[lv];
		}
		return sum;
	}
	template <typename DC_T>
	void output_stat(DC_T& dc){
		double isum = get_pr_sum();
		long nC = graph.get_num_combiners();
		long nS = graph.get_num_scatters();
		long nE = graph.get_num_edges();
		long allNC = 0;
		long allNS = 0;
		long allNE = 0;
		double sum = 0;
		const int root = 0;
		dc.refUnderlying().reduce(&isum, &sum, 1, MPI_DOUBLE, MPI_SUM, root);
		dc.refUnderlying().reduce(&nC, &allNC, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nS, &allNS, 1, MPI_LONG, MPI_SUM, root);
		dc.refUnderlying().reduce(&nE, &allNE, 1, MPI_LONG, MPI_SUM, root);
		if(dc.procID() == root){
			std::cout << std::setprecision(15) << "Sum of pr: "<< sum << std::endl;
			std::cout << "Number of combiners: "<< allNC << std::endl;
			std::cout << "Number of scatters: "<< allNS << std::endl;
			std::cout << "Number of agents: "<< allNC+allNS << std::endl;
			std::cout << "Number of edges: "<< allNE << std::endl;
			std::cout << "Value of equevalent cut: "<< (double)(allNC+allNS)/(double)(allNE) << std::endl;
		}
	}

private:
	vertex_data_t* pr;
	combine_data_t* pr_combine;
	scatter_data_t* pr_scatter;
}; 

#endif
