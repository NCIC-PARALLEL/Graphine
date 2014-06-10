#ifndef GRE_Util_FastForward_HPP
#define GRE_Util_FastForward_HPP
#include "Types.hpp"
#include "Config/SysConfig.hpp"
#include "Util/Locks.hpp"
#include <cstdlib>
#include <sched.h>
//
//Memory allocation should be further optimized!!!
//
namespace GRE{
	namespace Util{
		template <typename T>
		class FFQueue{
		public:
			typedef TicketLock lock_t;
			//typedef xcasLock lock_t;
		private:
			//Default construction is forbidden.
			FFQueue(){}
		public:
			FFQueue(T invalidV):invalidValue(invalidV),head(0), tail(0), size(0), sizeMask(0){}
			~FFQueue(){}
			void init(const size_t queueLen=1024, /*Set true only-if the invalidValue isn't 0*/const bool initValues=false){
				//1
				head = 0;
				headLock.init();
				tail = 0;
				tailLock.init();
				//2
				int log2Size = 1;
				while((0x1L << log2Size) < queueLen) log2Size++;
				size = 0x1L << log2Size;	
				sizeMask = size - 1;
				//queue = new T[size]invalidValue;
				//queue = (T*)malloc(sizeof(T)*size);
				queue = (T*)calloc(size, sizeof(T));
				if(initValues){for(int i=0; i<size; i++) queue[i]=invalidValue;}
			}
			void reset(){
				head = 0;
				tail = 0;
			}
			void destroy(){
				head = 0;
				tail = 0;
				size = 0;
				sizeMask = 0;
				free(queue);
				//delete[] queue;
			}
			inline bool isValid(T u){
				return (u != invalidValue);
			}
			//par
			inline T dequeue(){
				tailLock.spinLock();
				T tmp = queue[tail & sizeMask];
				if(isValid(tmp)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
				}
				tailLock.unLock();
				return tmp;
			}
			//par
			inline bool dequeue(T& value){
				tailLock.spinLock();
				value = queue[tail & sizeMask];
				if(isValid(value)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
					tailLock.unLock();
					return true;
				}
				tailLock.unLock();
				return false;
			}
			//par, maybe we should consider case "full".
			inline bool enqueue(const T& value){
				headLock.spinLock();
				if(!isValid(queue[head & sizeMask])){
					queue[head & sizeMask] = value;
					head++;
					headLock.unLock();
					return true;
				}
				headLock.unLock();
				return false;
			}
			//only single producer allowed
			inline T dequeue_serial(){
				T tmp = queue[tail & sizeMask];
				if(isValid(tmp)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
				}
				return tmp;
			}
			//only single producer allowed
			inline bool dequeue_serial(T& value){
				value = queue[tail & sizeMask];
				if(isValid(value)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
					return true;
				}
				return false;
			}
			//only single consumer allowed
			inline bool enqueue_serial(const T& value){
				if(!isValid(queue[head & sizeMask])){
					queue[head & sizeMask] = value;
					head++;
					return true;
				}
				return false;
			}
			//generally called by the producer!
			void waitForEmpty(){
				while(isValid(queue[(head-1) & sizeMask]))
					sched_yield();
			}
		private:
			//
			size_t head;
			lock_t headLock;
			byte* padding1[cacheLineSize-sizeof(size_t)-lock_t::size];
			size_t tail;
			lock_t tailLock;
			byte* padding2[cacheLineSize-sizeof(size_t)-lock_t::size];
			//
			size_t size;
			size_t sizeMask;
			T* queue;
			//
			const T invalidValue;
		};
/*
		template <typename T, T invalidValue>
		class FFQueue{
		public:
			typedef TicketLock lock_t;
		public:
			FFQueue():head(0), tail(0), size(0), sizeMask(0){}
			~FFQueue(){}
			void init(const size_t queueLen=1024){
				//1
				head = 0;
				headLock.init();
				tail = 0;
				tailLock.init();
				//2
				int log2Size = 1;
				while((0x1L << log2Size) < size) log2Size++;
				size = 0x1L << log2Size;	
				sizeMask = size - 1;
				queue = (T*)malloc(sizeof(T)*size);
			}
			void destroy(){
				head = 0;
				tail = 0;
				size = 0;
				sizeMask = 0;
				free(queue);
			}
			inline bool isValid(T& u){
				return (u != invalidValue);
			}
			//par
			inline T dequeue(){
				tailLock.spinLock();
				T tmp = queue[tail & sizeMask];
				if(isValid(tmp)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
				}
				tailLock.unLock();
				return tmp;
			}
			//par
			inline bool dequeue(T& value){
				tailLock.spinLock();
				value = queue[tail & sizeMask];
				if(isValid(value)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
					tailLock.unLock();
					return true;
				}
				tailLock.unLock();
				return false;
			}
			//par, maybe we should consider case "full".
			inline bool enqueue(const T& value){
				if(head - tail < size){
					headLock.spinLock();
					if(head - tail < size){
						queue[head & sizeMask] = value;
						head++;
						headLock.unLock();
						return true;
					}
					headLock.unLock();
				}
				return false;
			}
			//only single producer allowed
			inline T dequeue_serial(){
				T tmp = queue[tail & sizeMask];
				if(isValid(tmp)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
				}
				return tmp;
			}
			//only single producer allowed
			inline bool dequeue_serial(T& value){
				value = queue[tail & sizeMask];
				if(isValid(value)) {
					queue[tail & sizeMask] = invalidValue;
					tail++;
					return true;
				}
				return false;
			}
			//only single consumer allowed
			inline bool enqueue_serial(T& value){
				if(head - tail < size) {
					queue[head & sizeMask] = value;
					head++;
					return true;
				}
				return false;
			}
		private:
			//
			size_t head;
			lock_t headLock;
			byte* padding1[cacheLineSize-sizeof(size_t)-lock_t::size];
			size_t tail;
			lock_t tailLock;
			byte* padding2[cacheLineSize-sizeof(size_t)-lock_t::size];
			//
			size_t size;
			size_t sizeMask;
			T* queue;
		};
*/
	}//end Util
}
#endif
