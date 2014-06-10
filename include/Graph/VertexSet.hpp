#ifndef GRE_Graph_VtxSet_HPP
#define GRE_Graph_VtxSet_HPP

#include "Types.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Util/Bitmap.hpp"
#include "Util/Atomic.hpp"
#include "Util/Util.hpp"
#include <vector>

namespace GRE{
	namespace Graph{
		//Reenterent: support parallel operations
		//2*Phase:
		//1st:Current mode, i.e. Read-only
		//2nd:Next mode, i.e. Write-only
		class VertexSet{
		private:
			typedef GRE::Util::Bitmap Bitmap_t;
		private:
			bool isInitialized;

			//The active vertices can be either in form of bitmap or vector
			Bitmap_t vertexBitmap;

			int64_t capacity_;
			int64_t begin, end;
			int64_t end1;
			std::vector<lvid_t> vertexSet;
		public:
			VertexSet()
				:begin(0), end(0), end1(0), capacity_(0),vertexBitmap(), vertexSet(), isInitialized(false){}

			template <typename DistributedGraph_t>
			VertexSet(DistributedGraph_t& graph, const bool use_set = true)
				:begin(0), end(0), end1(0), capacity_(0),vertexBitmap(), vertexSet(){
				capacity_ = GRE::Util::align(graph.get_num_masters())+max(graph.get_num_scatters(), graph.get_num_combiners());
				vertexBitmap.resize(capacity_);
				if(use_set)
					vertexSet.resize(capacity_, -1);
				isInitialized = true;
			}

			VertexSet(const int _capacity, const bool use_set = true)
				:begin(0), end(0), end1(0), vertexBitmap(), vertexSet()	{
				capacity_ = _capacity;
				vertexBitmap.resize(capacity_);
				if(use_set)
					vertexSet.resize(capacity_, -1);
				isInitialized = true;
			}
			~VertexSet(){
				if(isInitialized) 
					destroy();	
			}

			void initialize(const int _capacity, const bool use_set = true){
				if(!isInitialized){
					begin = 0;
					end = 0;
					end1 = 0;
					capacity_ = _capacity;
					vertexBitmap.resize(capacity_);
					if(use_set)
						vertexSet.resize(capacity_, -1);
					isInitialized = true;
				}
			}
			template <typename DistributedGraph_t>
			void initialize(DistributedGraph_t& graph, const bool use_set = true){
				if(!isInitialized){
					begin = 0;
					end = 0;
					end1 = 0;
					//TODO:not safe in InOut edge mode.
					capacity_ = GRE::Util::align(graph.get_num_masters())+max(graph.get_num_scatters(), graph.get_num_combiners());
					vertexBitmap.resize(capacity_);
					if(use_set)
						vertexSet.resize(capacity_, -1);
					isInitialized = true;
				}
			}

			void destroy(){
				vertexBitmap.destroy();
				vertexSet.clear();
				isInitialized = false;
			}

			//Must be in next mode;
			inline void set(const lvid_t lvid)
			{
				vertexBitmap.set(lvid);
			}
			inline void unset(const lvid_t lvid)
			{
				vertexBitmap.unset(lvid);
			}
			inline bool test(const lvid_t lvid)
			{
				return vertexBitmap.test(lvid);
			}
			inline bool testAndSet(const lvid_t lvid)
			{
				return vertexBitmap.testAndSet(lvid);
			}


			//Must be in next mode;
			void insert(const lvid_t* buf, const int count)
			{
				//count should be >0, but not necessary.
				const int64_t baseIndex = AtomicFetchAdd(&end, count);
				for(int i = 0; i < count; i++){
					vertexSet[baseIndex+i] = buf[i];
				}
			}
			void insert(std::vector<lvid_t>& buf)
			{
				const int count = buf.size();
				//count should be >0, but not necessary.
				const int64_t baseIndex = AtomicFetchAdd(&end, count);
				for(int i = 0; i < count; i++){
					vertexSet[baseIndex+i] = buf[i];
				}
			}
			inline void insert(const lvid_t v)
			{
				const int64_t baseIndex = AtomicFetchAdd(&end, 1);
				vertexSet[baseIndex] = v;
			}

			//Must be in next mode;
			inline void setAndInsert(const lvid_t lvid)
			{
				if(!vertexBitmap.testAndSet(lvid))
				{
					const int64_t baseIndex = AtomicFetchAdd(&end, 1);
					vertexSet[baseIndex] = lvid;
				}
			}
			//It is better to implement an external buffer and use setAndInsert+insert.
			void setAndInsert(lvid_t* buf, const int count)
			{
				for(int i = 0; i < count; i++){
					const lvid_t lvid = buf[i];
					if(!vertexBitmap.testAndSet(lvid))
					{
						const int64_t baseIndex = AtomicFetchAdd(&end, 1);
						vertexSet[baseIndex] = lvid;
					}
				}
			}
			void setAndInsert(std::vector<lvid_t>& buf)
			{
				const int count = buf.size();	
				for(int i = 0; i < count; i++){
					const lvid_t lvid = buf[i];
					if(!vertexBitmap.testAndSet(lvid))
					{
						const int64_t baseIndex = AtomicFetchAdd(&end, 1);
						vertexSet[baseIndex] = lvid;
					}
				}
			}

			//Must be in current mode;
			int fetch(lvid_t* buf, const int count)
			{
				const int64_t baseIndex = AtomicFetchAdd(&begin, count);
				const int64_t left = end1 - baseIndex;
				if(left>0){
					const int64_t n = left > count? count:left;
					for(int i = 0; i < n; i++){
						buf[i] = vertexSet[baseIndex+i];
					}
					return n;
				} else
					return 0;
			}
			void refresh()
			{
				begin = 0;
				end1 = end;
				sync();
			}
			void clear(){
				begin = 0;
				end = 0;
				end1 = 0;
				vertexBitmap.clear();
				sync();
			}	
			void clear_set(){
				begin = 0;
				end = 0;
				end1 = 0;
				sync();
			}
			void clear_bitmap(){
				vertexBitmap.clear();
			}
			void reset(){
				clear();
			}

			inline int64_t size(){
				return (end);
			}
			inline int64_t num_ready(){
				return (end1);
			}
			inline int64_t num_new(){
				return (end - end1);
			}
			inline int64_t num_left(){
				return (end1-begin);
			}
			inline bool empty(){
				return (begin >= end1);
			}

			int64_t capacity(){
				return capacity_;
			}
			//
			//Fast transformation between set and bitmap
			//We introduce parallelism here.
			//
			//Must be in current mode;
			inline int fetchFromWord(lvid_t* buf, const lvid_t wordIdx){
				int count = 0;
				unsigned long tmpWord = vertexBitmap.getWord(wordIdx);
				lvid_t lvid = 64 * wordIdx;
				while(tmpWord){
					if(tmpWord & 0x1L){
						buf[count++] = lvid;	
					}
					lvid++;
					tmpWord>>=1;
				}
				return count;
			}

			void synchronizeBitmap2Set(){
				const int n = vertexBitmap.nWords();
				#pragma omp parallel
				{
					lvid_t buf[512];
					int count = 0;
					#pragma omp for
					for(int64_t i = 0; i<n; i++){
						unsigned long tmpWord = vertexBitmap.getWord(i);
						lvid_t lvid = 64 * i;
						while(tmpWord){
							if(tmpWord & 0x1L){
								buf[count++] = lvid;	
								if(count == 512){
									insert(buf, 512);
									count = 0;
								}
							}
							lvid++;
							tmpWord>>=1;
						}
					}
					insert(buf, count);
				}
			}

			//
			//Careful to use
			Bitmap_t& exportBitset(){
				return vertexBitmap;
			}
			std::vector<lvid_t>& exportSet(){
				return vertexSet;
			}
		private:
			template <typename T>
			inline T max(const T& t1, const T& t2){
				return (t1>t2? t1:t2);
			}
		};
		//
		struct SourceVertexSet{
			bool isComplete;
			std::vector<lvid_t> vertices;
		public:
			SourceVertexSet()
				:isComplete(false), vertices(){}

			bool isNULL(){
				if(isComplete || vertices.size()>0)
					return false;
				else
					return true;
			}
			void setNULL(){
				isComplete = false;
				vertices.clear();
			}

			bool isALL(){
				return isComplete;
			}
			void setALL(){
				isComplete = true;
			}

			void clear(){
				isComplete = false;
				vertices.clear();
			}
			void reset(){
				clear();
			}
		public:
			template <typename DistributedGraph>
			void loadFromFile(const std::string& prefix, DistributedGraph& graph){
				std::string filename = prefix + ".roots";
				std::ifstream infile(filename.c_str());
				int64_t n;
				infile >> n;
				for(int64_t i=0; i<n; i++){
					vid_t v;
					infile >> v;
					if(graph.is_local(v))
						vertices.push_back(graph.get_lvid_master(v));
				}
			}
			template <typename DistributedGraph>
			void addVertex(vid_t v, DistributedGraph& graph){
				if(graph.is_local(v))
					vertices.push_back(graph.get_lvid_master(v));
			}
			template <typename DistributedGraph>
			void addVertices(std::vector<vid_t>& src, DistributedGraph& graph, const size_t n=0){
				size_t cc;
				if(n==0 || n>src.size())
					cc = src.size();
				for(size_t i=0; i<cc; i++){
					if(graph.is_local(src[i]))
						vertices.push_back(graph.get_lvid_master(src[i]));
				}
			}
		};
	}//end Namespace-Graph
}
#endif
