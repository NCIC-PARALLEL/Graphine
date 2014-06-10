#ifndef GRE_Graph_DistributedGraph_HPP
#define GRE_Graph_DistributedGraph_HPP

#include <string>
#include "Types.hpp"
#include "Graph/LocalGraph.hpp"
#include "Graph/GraphIngress.hpp"
#include "Util/HashIndex.hpp"
#include "Util/Bitmap.hpp"
#include "Util/Util.hpp"
#include "Comm/BufferingComm.hpp"
#include "Comm/CommUnderlying.hpp"
#include <vector>
namespace GRE {
	namespace Graph {
		//1st Class to users.
		//Distributed Graph
		//
		class IngressParameters;
		template<typename T> class GraphIngress;

		template <typename LocalGraph>
		class DistributedGraph{
		public:
			typedef LocalGraph LocalGraph_t;
			typedef DistributedGraph<LocalGraph> SELF_t;
			friend class GraphIngress<SELF_t>;
		//public://types
			//typedef LocalGraph LocalGraph_t;
			//typedef VtxData VertexData_t;
			//typedef EdgeData EdgeData_t;
		public://interface
			//1:
			DistributedGraph()
				:isReady(false), subGraphID(0), graph(),vtxIndex(-1), scatterMaps(), combinerMaps(), rscatterMaps(), rcombinerMaps(){}
			~DistributedGraph(){}
			//2:Translate Rountines
			//2.1:global vertex id to local vertex id
			inline lvid_t get_lvid_master(const vid_t vid) const {
				return vid/gStat_nSubgraphs;
			}
			inline lvid_t get_lvid_scatter(const vid_t vid) {
				return vtxIndex.value(vid);
			}
			inline lvid_t get_lvid_combiner(const vid_t vid) {
				return vtxIndex.value(vid);
			}
			inline lvid_t get_lvid(const vid_t vid) const {
				if(vid%gStat_nSubgraphs == subGraphID)
					return vid/gStat_nSubgraphs;
				else
					return vtxIndex.value(vid);
			}

			//2.2:local vertex id to global vertex id
			//|...master...|...bi-agents...|...combiner...|...scatter...|
			inline vid_t get_vid(const lvid_t lvid) const {
				return lvid2vid[lvid];
			}
			inline int get_owner(vid_t v) const {
				return (v%gStat_nSubgraphs);
			}
			inline int is_local(const vid_t vid) const {
				return (get_owner(vid) == subGraphID);
			}
		
			inline size_t get_begin_lvid_master(){
				return (0);
			}
			inline size_t get_begin_lvid_agent(){
				return GRE::Util::align(lStat_nMasters);
			}
			inline size_t get_begin_lvid_scatter(){
				return GRE::Util::align(lStat_nMasters);
			}
			inline size_t get_begin_lvid_combiner(){
				return GRE::Util::align(lStat_nMasters);
			}

			inline size_t get_end_lvid_master(){
				return (lStat_nMasters);
			}
			inline size_t get_end_lvid_agent(){
				return (GRE::Util::align(lStat_nMasters)+lStat_nAgents);
			}
			inline size_t get_end_lvid_scatter(){
				return (GRE::Util::align(lStat_nMasters)+lStat_nScatters);
			}
			inline size_t get_end_lvid_combiner(){
				return (GRE::Util::align(lStat_nMasters)+lStat_nCombiners);
			}

			inline size_t get_num_masters(){
				return (lStat_nMasters);
			}
			inline size_t get_num_agents(){
				return (lStat_nAgents);
			}
			inline size_t get_num_scatters(){
				return (lStat_nScatters);
			}
			inline size_t get_num_combiners(){
				return (lStat_nCombiners);
			}
			inline size_t get_num_edges(){
				return lStat_nEdges;
			}
			//TODO:debug-only
			inline size_t get_num_vtx(){
				//return lStat_nMasters+GRE::Util::max(lStat_nScatters, lStat_nCombiners);
				return lStat_nMasters+lStat_nAgents;
			}

			//3:I/O
			void load(std::string prefix){
				//1: graph
				//
				//2: index-lvid2vid 
				//
				//3: index-vid2lvid
				//
				//4: scatter
				//
			}
			void save(std::string prefix){
				//1: graph
				//
				//2: index-lvid2vid
				//
				//3: index-vid2lvid
				//
				//4: scatter
				//
			}
			//interface
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, const std::string prefix){
				return ingress<Out>(comm, bufPool, prefix);
			}
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, const std::string prefix, 
				const std::string format, const std::string method){
				return ingress<Out>(comm, bufPool, prefix, format, method);
			}
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, GRE::Graph::IngressParameters& params){
				return ingress<Out>(comm, bufPool, params);
			}

			template<Dir dir>
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, const std::string prefix){
				IngressParameters params(prefix);
				return ingress<dir>(comm, bufPool, params);
			}
			template<Dir dir>
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, const std::string prefix, 
				const std::string format, const std::string method){
				IngressParameters params(prefix, format, method);
				return ingress<dir>(comm, bufPool, params);
			}

			//implementation
			template <Dir dir>
			bool ingress(GRE::COMM::Underlying& comm, GRE::Util::BufferPool& bufPool, GRE::Graph::IngressParameters& params){
				//set method
				GraphIngress<SELF_t>* ingress_p = new GraphIngress<SELF_t>(*this, comm, bufPool);
				//execute
				if(ingress_p->ingress<dir>(params)){
					isReady = true;
				}
				ingress_p->finalize();
				lStat_nEdges = graph.get_num_edges();
				subGraphID = comm.rank();
				return true;
			}
			//4: export
			inline LocalGraph_t& refLocalGraph(){
				return this->graph;
			}
			//5: property
			inline size_t getOutdegree(lvid_t lv){
				return outdegrees[lv];
			}
			template <typename T>
			void loadVertexProperty(T*& vData, const std::string filename){
				const size_t n=lStat_nMasters;
				vData = new T[n];
				std::ifstream infile(filename.c_str());
				if(!infile){
					std::cerr<<"Err in reading vertex peroperty from file"<<filename<<std::endl;
					exit(-1);
				}
				for(size_t v = 0; v < n; v++){
					infile >> vData[v];//TODO: catch exeptional case!
				}
			}
			template <typename T>
			void loadEdgeProperty(T*& eData, const std::string filename){
				const size_t n=lStat_nEdges;
				eData = new T[n];
				std::ifstream infile(filename.c_str());
				if(!infile){
					std::cerr<<"Err in reading edge peroperty from file"<<filename<<std::endl;
					exit(-1);
				}
				for(size_t e = 0; e < n; e++){
					infile >> eData[e];//TODO: catch exeptional case!
				}
			}
			//TODO:only for test
			template<typename T>
			void loadVertexProperty(T*& vData){
				const size_t n=lStat_nMasters;
				vData = new T[n];
				for(size_t v = 0; v < n; v++){
					vData[v] = rand()%65535;
					//vData[v] = 1;
				}
			}
			//TODO:only for test
			template<typename T>
			void loadEdgeProperty(T*& eData){
				const size_t n=lStat_nEdges;
				eData = new T[n];
				for(size_t e = 0; e < n; e++){
					eData[e] = rand()%65535;
					//eData[e] = 1;
				}
			}
			template<typename T>
			void unloadVertexProperty(T*& vData){
				delete[] vData;
			}
			template<typename T>
			void unloadEdgeProperty(T*& eData){
				delete[] eData;
			}

		private:
			typedef GRE::Util::FixedBitmap<MAX_NUM_MACHINES> FixedBitmap_t;
			typedef GRE::Util::Bitmap Bitmap_t;
			typedef GRE::Util::HashIndex<vid_t, lvid_t> HashIndex_t;
		private://data
			//Ready flag
			bool isReady;
			int subGraphID;
			//Global stats
			size_t gStat_nSubgraphs;
			size_t gStat_nVertices;
			size_t gStat_nEdges;
			//Local Stats
			size_t lStat_nScatters;
			size_t lStat_nCombiners;
			size_t lStat_nMasters;
			size_t lStat_nAgents;
			size_t lStat_nEdges;
			//
			//Local graph (Mandatory)
			LocalGraph_t graph;
			//vid <--> lvid index
			std::vector<vid_t> lvid2vid;
			//HashIndex_t vid2lvid;
			HashIndex_t vtxIndex;
			//HashIndex_t combinerIndex;
			//HashIndex_t MasterIndex;
			std::vector<int> outdegrees;
		public:
			//scatter map
			std::vector<Bitmap_t> scatterMaps;	
			std::vector<Bitmap_t> rscatterMaps;	
			std::vector<Bitmap_t> combinerMaps;	
			std::vector<Bitmap_t> rcombinerMaps;	
			//Local graph data (Optional)
			//std::vector<VertexData_t> vtxData;
			//std::vector<EdgeData_t> edgeData;
		};//end DG
	}//end Namespace-Graph
}
#endif
