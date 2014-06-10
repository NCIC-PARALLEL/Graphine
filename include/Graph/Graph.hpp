#ifndef GRE_Graph_HPP
#define GRE_Graph_HPP
#include "GraphTypesDef.hpp"
#include "UsrConfig.hpp"
namespace GRE{
	template <typename LocalGraph>
	class DistributedGraph
	{
	private:
		//global information
		size_t nVertices;
		size_t nEdges;
		int nSubgraphs;//should less than 64
		bool isLoaded;
		//local subgraph
		LocalGraph& localGraph;
	public:
		DistributedGraph():nVertices(0), nEdges(0), isLoaded(false){}
		DistributedGraph(LocalGraph& _localGraph):isLoaded(true), localGraph(_localGraph){
			//
		}
		~DistributedGraph(){}
	public:
		bool load(std::string graphName){
			//
		}
	};

	class LocalGraph
	{
	public:
	friend class GraphProperty;
	friend class GraphState;
	friend class DistributedGraph;
	friend class GraphImport;
	public:
	//types
	private:
		typedef GRE::Util::FixedBitmap<MAX_NUM_MACHINES> MachineMap_t;
		//ready flag
		bool isReady;
		//meta-data
		size_t nVertices;
		size_t nEdges;
		size_t nScatters;
		size_t nCombiners;
		//vid <--> lvid index
		std::vector<vid_t> lvid2vid;
		Util::HashIndex<vid_t, lvid_t> vid2lvid;
		//graph structure:Use CSR format
		std::vector<lidx_t> edgeIndex;
		std::vector<lvid_t> edgeList;
		//scatter agent distribution; 
		std::vector<MachineMap_t> scatterDistMap;
	public:
		localGraph(){}
		~localGraph(){}
	public:
		bool load(std::string prefix)
		{
			//load meta
			//
		}
		bool save(std::string prefix)
		{
			//
		}
	public:
		inline bool isScatter(lvid_t v)
		{
			return isScatter(lvid2vid[v]);	
		}
		inline bool isCombiner(lvid_t v)
		{
			return isCombiner(lvid2vid[v]);
		}
		inline bool isMaster(lvid_t v)
		{
			return isMaster(lvid2vid[v]);
		}
		inline machineMap_t getScatterList(lvid_t v)
		{
			return scatterDistMap[v];
		}
	public:
		//edgeIterator
		typedef lvid* edgeIterator;
		inline edgeIterator beginEdge(const lvid_t lvid){
			return &edgeList[lvid];
		}
		inline edgeIterator endEdge(const lvid_t lvid){
			return &edgeList[lvid];
		}
		//vtxIterator
		typedef lidx* vertexIterator;
		inline vertexIterator vertex_begin(){
			return &edge
		}
	};

	template <typename T>
	struct Edge{
		vid_t src;
		vid_t dest;
		T* data;	
	public:
		inline vid_t source()
		{
			return src;
		}
		inline vid_t target()
		{
			return dest;
		}
	};
	struct EdgeList{
		LocalGraph& localGraph;
		std::vector<Edge> edges;
	};

	template <typename T>
	struct Vertex{
		vid_t vid;
		T* data;
	public:
		inline vid_t id() const
		{
			return vid;
		}
		inline T* data()
		{
			return this->data;
		}
	};



	template <typename state_t>
	struct graphState{
		LocalGraph& graph;
		state_t*	states;
	};

	//Should not be changed.
	template <typename vertexProperty_t, typename edgeProperty_t>
	struct graphProperty{
		LocalGraph& graph;
		VertexProperty_t* vtxProperty;
		EdgeProperty_t*	edgeProperty;
	};
}
#endif
