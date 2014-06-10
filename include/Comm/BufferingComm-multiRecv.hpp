#ifndef GRE_COMM_BufComm_HPP
#define GRE_COMM_BufComm_HPP
#include "Types.hpp"
#include "Comm/CommUnderlying.hpp"
#include "Comm/CommFormat.hpp"
#include "Util/FIFOArray.hpp"
#include "Util/BufferPool.hpp"
#include "Util/Buffer.hpp"
#include <vector>
#include <assert.h>
namespace GRE{
	namespace COMM{
		class BufComm{
		private:
			typedef std::pair<byte*, size_t> packet_t;
		public:
			typedef GRE::Util::FIFOArray<packet_t,1024> outgoingQueue_t;
			typedef GRE::Util::FIFOArray<byte*,4096> incomingQueue_t;
			typedef GRE::Util::BufferPool BufferPool_t;
			typedef Underlying CommUnderlying_t;
		public:
			BufComm(CommUnderlying_t& _underlying, BufferPool_t& _pool, const bool do_initialize = false)
				:underlying(_underlying),pool(_pool),nActiveSending(0), nActiveRecving(0){
					if(do_initialize){
						initialize();
					}
				}
			~BufComm(){
				finalize();
			}
		private:
			//Part1: Buffer Infrastructure
			std::vector<outgoingQueue_t> outgoingQueues;
			incomingQueue_t incomingQueue;

			std::vector<byte*> outgoingBuf;
			std::vector<byte*> incomingBuf;

			std::vector<int> outgoingTag;
			std::vector<int> incomingTag;

			//Part2: Data Transfering
			CommUnderlying_t& underlying;	
			//Part3:
			BufferPool_t& pool;
			//Part4
			int nActiveSending;
			int nActiveRecving;
			//Part5: should not go here!
			std::vector<int> backup;
		public://Interface
			void initialize(){
				outgoingQueues.resize(underlying.size());
				outgoingBuf.resize(underlying.size(),NULL);
				outgoingTag.resize(underlying.size(),0);
				incomingBuf.resize(underlying.size(),NULL);
				incomingTag.resize(underlying.size(),0);
			}
			void finalize(){
				outgoingBuf.clear();
				outgoingTag.clear();
				outgoingQueues.clear();
				incomingBuf.clear();
				incomingTag.clear();
				backup.clear();
			}
			template <typename Buffer_t>
			inline void isend(Buffer_t& buf,const int dest)
			{
				packet_t pkt;
				pkt.first = buf.address();
				pkt.second = buf.length();
				//byte* addr = buf.address();
				//size_t len = buf.length();
				if(outgoingBuf[dest]){//not NULL
					if(!underlying.testSend(dest)){//busy!
						//outgoingQueues[dest].put(make_pair(addr, len));
						outgoingQueues[dest].put(pkt);
						return;
					}
					pool.put_address(outgoingBuf[dest]);
				} else {
					nActiveSending++;
				}
			
				if(outgoingQueues[dest].empty()){//bypass and send imediately
					outgoingBuf[dest] = pkt.first;
					isend_internal(pkt.first, pkt.second, dest, outgoingTag[dest]++);
					//outgoingBuf[dest] = addr;
					//isend_internal(addr, len, dest);
				} else {
					const packet_t pkt1 = outgoingQueues[dest].get();
					//outgoingQueues[dest].put(make_pair(addr, len));
					outgoingQueues[dest].put(pkt);
					isend_internal(pkt1.first, pkt1.second,dest, outgoingTag[dest]++);
					outgoingBuf[dest] = pkt1.first;
				}
			}

			bool checkAllSend()
			{
				if(nActiveSending == 0) return false; //no sending active.
				int nActive = 0;
				for(int i=0; i<underlying.size(); i++){
					if(outgoingBuf[i]){// !=NULL
						if(underlying.testSend(i)){//finished
							pool.put_address(outgoingBuf[i]);
						} else {
							nActive++;
							continue;
						}
						if(!outgoingQueues[i].empty()){
							const packet_t pkt = outgoingQueues[i].get();
							isend_internal(pkt.first, pkt.second, i, outgoingTag[i]++);
							outgoingBuf[i] = pkt.first;
							nActive++;
						} else {
							outgoingBuf[i] = NULL;
						}
					}
				}
				nActiveSending = nActive;
				return (nActive != 0);
			}

			bool checkSend(int dest){
				if(!outgoingBuf[dest]) return false;//buf is NULL, i.e. not active
				if(underlying.testSend(dest)){//finished
					pool.put_address(outgoingBuf[dest]);
				} else 
					return true;
				if(!outgoingQueues[dest].empty()){
					const packet_t pkt = outgoingQueues[dest].get();
					isend_internal(pkt.first, pkt.second, dest, outgoingTag[dest]++);
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
				if(incomingBuf[src]){
					if(!underlying.testRecv(src)){//busy
						return true;
					}
					incomingQueue.put(incomingBuf[src]);
				} else {
					nActiveRecving++;
				}
				//lauch a new irecv from src
				incomingBuf[src] = pool.get_address();
				assert(incomingBuf[src]!=NULL);
				if(!incomingBuf[src]){//Err: no buffer to receive data...
					return false;
				}
				return underlying.irecv(incomingBuf[src], pool.bufSize(), src, incomingTag[src]++);
			}

			//Return true means there are still jobs on.
			bool checkAllRecv(const bool restart=true){
				if(nActiveRecving==0) return false;//
				int nActive = 0;
				for(int i = 0; i < underlying.size(); i++){
					if(incomingBuf[i]){
						if(underlying.testRecv(i)){//finished
							incomingQueue.put(incomingBuf[i]);
							if(restart){
								incomingBuf[i] = pool.get_address();
								assert(incomingBuf[i]!=NULL);
								underlying.irecv(incomingBuf[i], pool.bufSize(), i, incomingTag[i]++);//restart
								nActive++;
							} else
								incomingBuf[i] = NULL;
						} else
							nActive++;
					}
				}
				nActiveRecving = nActive;
				return (nActive != 0);
			}

			//Return true means there are still jobs on.
			bool checkRecv(int src, const bool restart = true){
				if(!incomingBuf[src]) return false;
				if(underlying.testRecv(src)){//finished
					incomingQueue.put(incomingBuf[src]);
					if(restart){
						incomingBuf[src] = pool.get_address();
						assert(incomingBuf[src]!=NULL);
						underlying.irecv(incomingBuf[src], pool.bufSize(), src, incomingTag[src]++);//restart
						return true;
					}
					incomingBuf[src] = NULL;
					nActiveRecving--;
					return false;
				}
				return true;
			}
			inline void cancelRecv(int src)
			{
				if(incomingBuf[src]){
					if(!underlying.testRecv(src))
						underlying.cancelRecv(src);
					pool.put_address(incomingBuf[src]);
					incomingTag[src] = 0;
					incomingBuf[src] = NULL;
					nActiveRecving--;
				}
			}
			inline void cancelSend(int dest)
			{
				if(outgoingBuf[dest]){
					if(!underlying.testSend(dest))
						underlying.cancelSend(dest);
					pool.put_address(outgoingBuf[dest]);
					outgoingBuf[dest] = NULL;
					outgoingTag[dest] = 0;
					nActiveSending--;
				}
			}
			void cancelAllRecv()
			{
				for(int i=0; i<underlying.size(); i++){
					if(incomingBuf[i]){
						if(!underlying.testRecv(i))
							underlying.cancelRecv(i);
						pool.put_address(incomingBuf[i]);
						incomingTag[i] = 0;
						incomingBuf[i] = NULL;
						nActiveRecving--;
					}
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
						outgoingTag[i] = 0;
					}
				}
				nActiveSending = 0;
			}
			bool launchAllRecv()
			{
				for(int i=0; i<underlying.size(); i++){
					if(i != underlying.rank())//don't receive from myself!
						if(!irecv(i)) return false;//Error!
				}
				return true;
			}
			int nActiveRecv(){
				return nActiveRecving;
			}
			int nActiveSend(){
				return nActiveSending;
			}
			//We must make the incoming queue visible to outside
			incomingQueue_t& getIncomingQueue(){
				return incomingQueue;
			}
			void exportIncomingQueue(incomingQueue_t* queue){
				queue = &incomingQueue;
			}
		private:
			inline void isend_internal(byte* addr, const size_t count, const int dest, const int tag = 0){
				underlying.isend(addr, count, dest, tag);
			}
		};
	}//end COMM
}
#endif
