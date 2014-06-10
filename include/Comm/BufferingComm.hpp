#ifndef GRE_COMM_BufComm_HPP
#define GRE_COMM_BufComm_HPP
#include "Types.hpp"
#include "Comm/CommUnderlying.hpp"
#include "Comm/Protocol.hpp"
#include "Util/FIFOArray.hpp"
#include "Util/BufferPool.hpp"
#include "Util/Buffer.hpp"
#include <vector>
#include <assert.h>
namespace GRE{
	namespace COMM{
		template<typename BufferPool_t>
		class BufComm{
		private:
			typedef std::pair<byte*, size_t> packet_t;
		public:
			typedef GRE::Util::FIFOArray<packet_t,2048> outgoingQueue_t;
			typedef GRE::Util::FIFOArray<byte*,8192> incomingQueue_t;
			typedef Underlying CommUnderlying_t;
		public:
			BufComm(CommUnderlying_t& _underlying, BufferPool_t& _pool, const bool do_initialize = true)
				:underlying(_underlying),pool(_pool),nActiveSending(0), nActiveRecving(0){
					if(do_initialize){
						initialize();
					} else {
						initialized = false;
					}
				}
			~BufComm(){
				finalize();
			}
		private:
			bool initialized;
			//Part1: Buffer Infrastructure
			std::vector<outgoingQueue_t> outgoingQueues;
			incomingQueue_t incomingQueue;

			std::vector<byte*> outgoingBuf;
			std::vector<bool> incomingFlags;
			byte* incomingBuf;

			//Part2: Data Transfering
			CommUnderlying_t& underlying;	
			//Part3:
			BufferPool_t& pool;
			//Part4
			int nActiveSending;
			int nActiveRecving;
			//Part5: should not go here!
		public://Interface
			void initialize(){
				if(!initialized){
					outgoingQueues.resize(underlying.size());
					outgoingBuf.resize(underlying.size(),NULL);
					incomingFlags.resize(underlying.size(),false);
					nActiveSending = 0;
					nActiveRecving = 0;
					incomingBuf = NULL;
				}
			}
			void reset(){
				if(initialized){
					if(outgoingQueues.size() < underlying.size())
						outgoingQueues.resize(underlying.size());
					for(int i=0; i<underlying.size(); i++){
						outgoingQueues[i].reset();
					}
					if(outgoingBuf.size() < underlying.size()){
						outgoingBuf.resize(underlying.size(),NULL);
						incomingFlags.resize(underlying.size(), false);
					} else {
						for(int i=0; i<underlying.size(); i++){
							outgoingBuf[i] = NULL;
							incomingFlags[i] = false;
						}
					}
					incomingBuf = NULL;
					nActiveSending = 0;
					nActiveRecving = 0;
				}
			}
			void finalize(){
				#ifdef Debug
				for(int i=0; i<underlying.size(); i++){
					assert(!outgoingBuf[i] && outgoingQueues[i].empty());
				}
				assert(!incomingBuf && incomingQueue.empty());
				#endif
				outgoingBuf.clear();
				outgoingQueues.clear();
			}
			template <typename Buffer_t>
			inline void isend(Buffer_t& buf,const int dest)
			{
				byte* addr = buf.address();
				size_t len = buf.length();
				if(outgoingBuf[dest]){//not NULL
					if(!underlying.testSend(dest)){//busy!
						packet_t pkt(addr, len);
						outgoingQueues[dest].put(pkt);
						return;
					}
					pool.put_address(outgoingBuf[dest]);
				} else {
					nActiveSending++;
				}
			
				if(outgoingQueues[dest].empty()){//bypass and send imediately
					outgoingBuf[dest] = addr;
					isend_internal(addr, len, dest);
				} else {
					packet_t pkt(addr, len);
					outgoingQueues[dest].put(pkt);
					const packet_t pkt1 = outgoingQueues[dest].get();
					isend_internal(pkt1.first, pkt1.second, dest);
					outgoingBuf[dest] = pkt1.first;
				}
			}

			inline bool checkAllSend()
			{
				if(nActiveSending==0) return false; 
				int nActive = 0;
				for(int i=0; i<underlying.size(); i++){
					if(outgoingBuf[i]){// !=NULL
						assert(i != underlying.rank());
						if(underlying.testSend(i)){//finished
							pool.put_address(outgoingBuf[i]);
						} else {
							nActive++;
							continue;
						}
						if(!outgoingQueues[i].empty()){
							const packet_t pkt = outgoingQueues[i].get();
							isend_internal(pkt.first, pkt.second, i);
							outgoingBuf[i] = pkt.first;
							nActive++;
						} else {
							outgoingBuf[i] = NULL;
						}
					}
				}
				nActiveSending = nActive;
				return (bool)(nActive != 0);
			}

			bool checkSend(int dest){
				if(!outgoingBuf[dest]) return false;//buf is NULL, i.e. not active
				if(underlying.testSend(dest)){//finished
					pool.put_address(outgoingBuf[dest]);
				} else 
					return true;
				if(!outgoingQueues[dest].empty()){
					const packet_t pkt = outgoingQueues[dest].get();
					isend_internal(pkt.first, pkt.second, dest);
					outgoingBuf[dest] = pkt.first;
					return true;
				} else {
					outgoingBuf[dest] = NULL;//Done!
					nActiveSending--;
					return false;
				}
			}

			//Return true means recv is on.
			inline bool irecv(int src)
			{
				assert(src!=underlying.rank());//don's send to myself!
				if(!incomingFlags[src]){//
					nActiveRecving++;
					incomingFlags[src] = true;
				}
				//check and recv
				if(incomingBuf){
					if(!underlying.testRecv()){//busy
						return true;
					}
					incomingQueue.put(incomingBuf);
				} 
				//lauch a new irecv from src
				incomingBuf = pool.get_address();
				assert(incomingBuf!=NULL);
				if(!incomingBuf){//Err: no buffer to receive data...
					return false;
				}
				return underlying.irecv(incomingBuf, pool.bufSize());
			}

			bool checkAllRecv(){
				if(!incomingBuf) return false;
				if(underlying.testRecv()){//finished
					byte* addr=pool.get_address();
					if(addr){
						incomingQueue.put(incomingBuf);
						incomingBuf = addr;
						underlying.irecv(incomingBuf, pool.bufSize());//restart
					}
				}
				return true;
			}
			//Return true means there are still jobs on.
			bool checkAllRecv(bool& relaunched){
				relaunched = true;
				if(nActiveRecving==0)
					return false;
				if(incomingBuf){
					if(underlying.testRecv()){//finished
						incomingQueue.put(incomingBuf);
						incomingBuf = NULL;
						relaunched = resumeRecv_();
					}
				} else {
					relaunched = resumeRecv_();
				}
				return true;
			}
			bool resumeRecv(){
				if(incomingBuf) return true;//yet recving
				if(nActiveRecving==0) return false;//no need to launch
				return resumeRecv_();
			}
			//private:
			inline bool resumeRecv_(){
				byte* addr = pool.get_address();
				if(addr){
					incomingBuf = addr;
					underlying.irecv(incomingBuf, pool.bufSize());//restart
					return true;
				}
				return false;
			}
	
			bool checkRecv(int src){
				if(!incomingFlags[src]) return false;
				if(underlying.testRecv()){//finished
					incomingQueue.put(incomingBuf);
					incomingBuf = pool.get_address();
					assert(incomingBuf);
					underlying.irecv(incomingBuf, pool.bufSize());//restart
				}
				return true;
			}

			inline void cancelRecv(int src)
			{
				if(incomingFlags[src]){
					incomingFlags[src] = false;
					nActiveRecving--;
				}
				if(nActiveRecving==0){
					if(incomingBuf){
						if(!underlying.testRecv())
							underlying.cancelRecv();
						pool.put_address(incomingBuf);
						incomingBuf = NULL;
					}
				}
			}
			inline void cancelSend(int dest)
			{
				if(outgoingBuf[dest]){
					if(!underlying.testSend(dest))
						underlying.cancelSend(dest);
					pool.put_address(outgoingBuf[dest]);
					outgoingBuf[dest] = NULL;
					nActiveSending--;
				}
			}
			void cancelAllRecv()
			{
				for(int i=0; i<underlying.size(); i++){
					if(incomingFlags[i]){
						incomingFlags[i] = false;
					}
				}
				nActiveRecving = 0;
				if(incomingBuf){
					if(!underlying.testRecv())
						underlying.cancelRecv();
					pool.put_address(incomingBuf);
					incomingBuf = NULL;
				}
			}
			void cancelAllSend()
			{
				for(int i=0; i<underlying.size(); i++){
					if(outgoingBuf[i]){
						if(!underlying.testSend(i))
							underlying.cancelSend(i);
						pool.put_address(outgoingBuf[i]);
						outgoingBuf[i] = NULL;
					}
				}
				nActiveSending = 0;
			}
			bool launchAllRecv()
			{
				nActiveRecving = underlying.size()-1;
				if(nActiveRecving==0) return true;
				for(int i = 0; i<underlying.size(); i++){
					if(i!=underlying.rank())
						incomingFlags[i] = true;
				}
				//lauch a new irecv from src
				incomingBuf = pool.get_address();
				assert(incomingBuf);
				if(!incomingBuf){//Err: no buffer to receive data...
					return false;
				}
				return underlying.irecv(incomingBuf, pool.bufSize());
			}
			inline int nActiveRecv(){
				return nActiveRecving;
			}
			inline int nActiveSend(){
				return nActiveSending;
			}
			inline int nProcs(){
				return underlying.size();
			}
			inline int procID(){
				return underlying.rank();
			}
			inline void barrier(){
				underlying.barrier();
			}
			//We must make the incoming queue visible to outside
			incomingQueue_t& getIncomingQueue(){
				return incomingQueue;
			}
			//same with the above
			incomingQueue_t& getRefIncomingQueue(){
				return incomingQueue;
			}
			void exportIncomingQueue(incomingQueue_t* queue){
				queue = &incomingQueue;
			}
			CommUnderlying_t& getRefUnderlying(){
				return underlying;
			}
			CommUnderlying_t& refUnderlying(){
				return underlying;
			}
			int voteToHalt(const bool isActive=false){
				return underlying.voteToHalt(isActive);
			}
		private:
			inline void isend_internal(byte* addr, const size_t count, const int dest, const int tag = 0){
				underlying.isend(addr, count, dest, tag);
			}
		};
	}//end COMM
}
#endif
