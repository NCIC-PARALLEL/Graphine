#ifndef _GRE_VertexCode_HPP
#define _GRE_VertexCode_HPP
#include "Types.hpp"
#include "Message.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Util/Util.hpp"

namespace GRE{
	template<typename scatter_data_t_, typename combine_data_t_, typename vertex_data_t_, typename DistributedGraph_t>
	class VertexCode{
	public:
		typedef scatter_data_t_ scatter_data_t;
		typedef combine_data_t_ combine_data_t;
		typedef vertex_data_t_ vertex_data_t;

		typedef scatter_data_t ScatterData_t;
		typedef combine_data_t CombineData_t;
		typedef vertex_data_t VertexData_t;
	public:

		//////////////////////////////////////////////////////////////////////////////////////////
		//Vertex-Code Interface
		////
		//Init
		bool init(const lvid_t lv){
			std::cerr<<"Fetal Error! Init primitive is not implemented."<<std::endl;
			return false;
		}

		//Scatter
		template <typename Ctx_t, typename EdgeIterator_t>
		void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
			std::cerr<<"Fetal Error! Scatter primitive is not implemented."<<std::endl;
		}

		//Combine
		template <typename Lock_t>
		bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
			std::cerr<<"Fetal Error! Combine primitive is not implemented."<<std::endl;
			return false;
		}

		//Apply
		bool apply(const lvid_t lv){
			std::cerr<<"Fetal Error! Apply primitive is not implemented."<<std::endl;
			return false;
		}

		//Assert-to-Halt
		bool assert_to_halt(const lvid_t lv){return true;}

		/////////////////////////////////////////////
		//Extension-1: Aggregate
		//Spec: called by each active vertex
		//local operation
		void aggregate(lvid_t lv){}

		/////////////////////////////////////////////
		//Extension-2: reduce
		//Spec: called once (i.e., by all active vertices) per superstep
		//conceptually global operation across cluster
		template <typename DC_t>
		void reduce(DC_t& dc){}

		/////////////////////////////////////////////
		//Extension-2: reduce
		//Spec: called once (i.e., by all active vertices) per superstep
		//conceptually global operation across cluster
		template <typename DC_t>
		void output_statistics(DC_t& dc){}


		public:
		//Default implementation. Should not modify in derived classes.
		template <typename Ctx_t>
		void copy_scatter_data(Ctx_t& ctx, scatter_data_t& data, const lvid_t lsrc){
			data = scatter_data[lsrc];
		}
		inline void copy_scatter_data(scatter_data_t& data, const lvid_t lsrc){
			data = scatter_data[lsrc];
		}

		template <typename Ctx_t>
		void update_scatter(Ctx_t& ctx, const lvid_t ldst, const scatter_data_t& data){
			scatter_data[ldst]=data;
		}
		inline void update_scatter(const lvid_t ldst, const scatter_data_t& data){
			scatter_data[ldst]=data;
		}

		template <typename Ctx_t>
		void copy_combine_data(Ctx_t& ctx, combine_data_t& data, const lvid_t lv){
			data = combine_data[lv];
		}
		inline void copy_combine_data(combine_data_t& data, const lvid_t lv){
			data = combine_data[lv];
		}

		template <typename Ctx_t>
		void reset_combiner(Ctx_t& ctx, const lvid_t ldst){
			//must define!!!
			std::cerr<<"Fetal Error! reset combiner primitive is not implemented."<<std::endl;
		}
		inline void reset_combiner(const lvid_t ldst){
			//must define!!!
			std::cerr<<"Fetal Error! reset combiner primitive  is not implemented."<<std::endl;
		}

		template <typename Ctx_t>
		inline void reset_scatter(Ctx_t& ctx, const lvid_t lsrc){
			//in default, do nothing.
		}
		inline void reset_scatter(const lvid_t lsrc){
			//in default, do nothing.
		}

	protected:
		scatter_data_t* scatter_data;
		combine_data_t* combine_data;
		vertex_data_t* vertex_data;

		DistributedGraph_t& graph;
	protected:
		//forbid explicit assignment.
		VertexCode(DistributedGraph_t& _graph)
			:graph(_graph),scatter_data(NULL),combine_data(NULL),vertex_data(NULL)
		{}
	public:
		~VertexCode(){}

	public:
		void allocScatterData(scatter_data_t*& ptr){
			assert(!ptr);
			const long n = GRE::Util::align(graph.get_num_masters())+graph.get_num_scatters();
			ptr = new scatter_data_t[n];
			scatter_data = ptr;
		}
		void allocScatterData(scatter_data_t*& ptr, const scatter_data_t iv){
			assert(!ptr);
			const long n = GRE::Util::align(graph.get_num_masters())+graph.get_num_scatters();
			ptr = new scatter_data_t[n];
			scatter_data = ptr;
			#pragma omp parallel for
			for(long i=0; i<n; i++){
				scatter_data[i] = iv;
			}
		}
		void allocCombineData(combine_data_t*& ptr){
			assert(!ptr);
			const long n = GRE::Util::align(graph.get_num_masters())+graph.get_num_combiners();
			ptr = new combine_data_t[n];
			combine_data = ptr;
		}
		void allocCombineData(combine_data_t*& ptr, const combine_data_t iv){
			assert(!ptr);
			const long n = GRE::Util::align(graph.get_num_masters())+graph.get_num_combiners();
			ptr = new combine_data_t[n];
			combine_data = ptr;
			#pragma omp parallel for
			for(long i=0; i<n; i++){
				combine_data[i] = iv;
			}
		}
		void allocVertexData(vertex_data_t*& ptr){
			if(!ptr)//have not been allocated.
				ptr = new vertex_data_t[graph.get_num_masters()];
			vertex_data = ptr;
		}
		void allocVertexData(vertex_data_t*& ptr, const vertex_data_t iv){
			const long n=graph.get_num_masters();
			if(!ptr)//have not been allocated.
				ptr = new vertex_data_t[n];
			vertex_data = ptr;
			#pragma omp parallel for
			for(long i=0; i<n; i++){
				vertex_data[i] = iv;
			}
		}


		void setScatterData(const scatter_data_t* ptr){
			scatter_data = ptr;
		}
		void setCombineData(const combine_data_t* ptr){
			combine_data = ptr;
		}
		void setVertexData(const vertex_data_t* ptr){
			vertex_data = ptr;
		}

		void freeScatterData(scatter_data_t*& ptr){
			if(ptr){
				delete[] ptr;
				ptr = NULL;
			}
			scatter_data = NULL;
		}
		void freeCombineData(combine_data_t*& ptr){
			if(ptr){
				delete[] ptr;
				ptr = NULL;
			}
			combine_data = NULL;
		}
		void freeVertexData(vertex_data_t*& ptr){
			if(ptr){
				delete[] ptr;
				ptr = NULL;
			}
			vertex_data = NULL;
		}
	};
}


#endif
