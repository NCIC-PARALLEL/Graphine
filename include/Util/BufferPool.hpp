#ifndef GRE_COMM_BufferPool_HPP
#define GRE_COMM_BufferPool_HPP
#include "Types.hpp"
#include "Config/UsrConfig.hpp"
#include "Util/Locks.hpp"
#include <cstdlib>
#include <vector>
namespace GRE{
	namespace Util{
		class BufferPool{	
		public:
			BufferPool():initialized(false), isOwner(false), pool(NULL), _bufSize(0), _bufCount(0), _capacity(0){}
			virtual ~BufferPool(){
				if(initialized)
					destroy();
			}
			bool config(const size_t bufCount, const size_t bufSize = middleBufSize, const bool do_initialize=false){
				this->_bufCount = bufCount;
				this->_bufSize = bufSize;
				this->_capacity = bufCount * bufSize;
				if(do_initialize){
					pool = (byte*)malloc(_capacity*sizeof(byte));
					if(!pool) return false;// fail to allocate!
					freeList.resize(bufCount);
					for(size_t i=0; i<bufCount; i++){
						freeList[i] = i;
					}
					this->isOwner = true;
					this->initialized = true;
				}
				return true;
			}
			bool config()
			{
				return config(defaultBufPoolSize, middleBufSize, false);
			}
			bool resetBufSize(const size_t bufSize){
				if(!initialized) return false;
				if(this->_bufSize != bufSize){
					_bufSize = bufSize;
					_bufCount = _capacity/bufSize;
					freeList.clear();
					freeList.resize(_bufCount);
					for(size_t i=0; i<_bufCount; i++){
						freeList[i] = i;
					}
				}
				return true;
			}
			bool resetBufCount(const size_t bufCount){
				if(!initialized) return false;
				if(this->_bufCount != bufCount){
					_bufCount = bufCount;
					_bufSize = _capacity/bufCount;
					freeList.clear();
					freeList.resize(_bufCount);
					for(size_t i=0; i<_bufCount; i++){
						freeList[i] = i;
					}
				}
				return true;
			}

			bool initialize()
			{
				if(!initialized){
					pool = (byte*)malloc(_bufCount*_bufSize);
					if(!pool) return false;// fail to allocate!
					isOwner = true;
					freeList.resize(_bufCount);
					for(size_t i=0; i<_bufCount; i++){
						freeList[i] = i;
					}
					initialized = true;
				}
				return true;
			}
			inline bool isInitialized(){
				return initialized;
			}

			bool inherit(BufferPool& anotherPool, const size_t newBufCount=0){
				if(initialized) //Can't inherit if have been initialized.
					return false;
				pool = anotherPool.base();	
				isOwner = false;
				_capacity = anotherPool.capacity();
				_bufCount = (newBufCount==0? anotherPool.bufCount() : newBufCount);
				_bufSize = _capacity/_bufCount;
				freeList.clear();
				freeList.resize(_bufCount);
				for(size_t i=0; i<_bufCount; i++){
					freeList[i] = i;
				}
				initialized = true;
				return true;
			}

			void destroy()
			{
				if(initialized){
					initialized = false;
					//assert(_bufCount == freeList.size());
					_bufCount = 0;
					_bufSize = 0;
					_capacity = 0;
					freeList.clear();
					if(isOwner == true)
						free(pool);
				}
			}
			void finalize()
			{
				this->destroy();
			}
			inline size_t get(){
				if(freeList.empty())
					return -1;
				else{
					size_t tmp = freeList.back();
					freeList.pop_back();
					return tmp;
				}
			}
			inline byte* get_address()
			{
				if(freeList.empty())
					return NULL;
				else{
					const size_t tmp = freeList.back();
					freeList.pop_back();
					return (pool+_bufSize*tmp);
				}
			}
			inline void put(const size_t bufID){
				freeList.push_back(bufID);
			}
			inline void put_address(const byte* addr){
				freeList.push_back(bufID(addr));
			}

			void reset(){
				freeList.resize(_bufCount);
				for(size_t i=0; i<_bufCount; i++){
					freeList[i] = i;
				}
			}

			inline size_t bufSize(){
				return _bufSize;
			}
			inline size_t bufCount(){
				return _bufCount;
			}
			inline size_t freeCount(){
				return freeList.size();
			}
			inline size_t capacity(){
				return _capacity;
			}

			//
			byte* address(const size_t bufID){
				if(bufID < 0 || bufID >= _bufCount)
					return NULL;
				return (pool+_bufSize*bufID);
			}
			//
			size_t bufID(const byte* addr){
				if(addr < pool || addr >= pool+_capacity || (addr - pool)%_bufSize!=0)
					return -1;
				return (size_t)((addr - pool)/_bufSize);
			}
			//Don't use it explicitly!
			byte* base(){
				return pool;
			}
		private:
			bool initialized;
			bool isOwner;
			std::vector<size_t> freeList;
			byte* pool;

			//capacity = bufSize * (maxBufID+1)
			size_t _bufSize;
			size_t _bufCount;
			size_t _capacity;
		};
		//parBufferPool
		class parBufferPool:public BufferPool{
		private:
			typedef GRE::Util::TicketLock lock_t;
			lock_t lock;
		public:
			parBufferPool():lock(){}
			~parBufferPool(){}
			size_t get(){
				size_t ret;
				lock.spinLock();
				ret = BufferPool::get();
				lock.unLock();
				return ret;
			}
			byte* get_address(){
				byte* ret;
				lock.spinLock();
				ret = BufferPool::get_address();
				lock.unLock();
				return ret;
			}
			void put(const size_t bufID){
				lock.spinLock();
				BufferPool::put(bufID);
				lock.unLock();
			}
			void put_address(const byte* addr){
				lock.spinLock();
				BufferPool::put_address(addr);
				lock.unLock();
			}
		};
	}//end Comm
}
#endif
