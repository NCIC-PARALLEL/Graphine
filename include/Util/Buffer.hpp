#ifndef GRE_Util_Buffer_HPP
#define GRE_Util_Buffer_HPP
#include "Types.hpp"
#include "Util/Locks.hpp"
namespace GRE{
	namespace Util{
		//
		//Internal View:
		//||Header|data|data|data|data|...
		//
		template <typename Header_t, typename Data_t>
		class FormatedBuffer{
		public:
			typedef Data_t* data_iterator;
		public:
			FormatedBuffer():header(NULL),data(NULL),size_(0){}
			FormatedBuffer(byte* addr, size_t capacity, const bool do_format=false){
				header = (Header_t*)addr;
				data = (Data_t*)(addr + sizeof(Header_t));
				size_ = (capacity - sizeof(Header_t))/sizeof(Data_t);
				if(do_format)
					format();
			}
			~FormatedBuffer(){}

			void open(byte* addr, size_t capacity, const bool do_format=false){
				header = (Header_t*)addr;
				data = (Data_t*)(addr + sizeof(Header_t));
				size_ = (capacity - sizeof(Header_t))/sizeof(Data_t);
				if(do_format)
					format();
			}
			void close(){
				header = NULL;
				data = NULL;
				size_ = 0;
			}

			//obsolete, please use open and than format, or open with do_format.
			void format(const byte* addr, const size_t capacity){
				header = (Header_t*)addr;
				header->clear();
				data = (Data_t*)(addr + sizeof(Header_t));
				size_ = (capacity - sizeof(Header_t))/sizeof(Data_t);
			}

			void format(){
				header->clear();
			}

			inline bool ready(){
				return (header!=NULL);
			}
			inline bool empty(){
				return (header->count == 0);
			}
			inline bool full(){
				return (header->count == size_); 
			}
			inline size_t capacity(){
				return size_*sizeof(Data_t);
			}
			inline size_t size(){
				return size_;
			}
			inline void clear(){
				header->clear();	
			}

			inline size_t length(){
				return (sizeof(Header_t)+size_*sizeof(Data_t));
			}
			inline size_t count(){
				return header->count;
			}
		
			//put
			//Require: Data_t must be assignable!
			inline bool push_back(const Data_t& val){
				if(!this->full()){
					data[header->count++] = val;
					return true;// success
				} else {
					return false; // failure
				}
			}
			inline bool pop_back(Data_t& val){
				if(!this->empty()){
					val = data[header->count-1];
					header->count--;
					return true;// success
				} else {
					return false; // failure
				}
			}
			inline void pop_back(){
				if(!this->empty()){
					header->count--;
				}
			}
			inline bool back(Data_t& val){
				if(!this->empty()){
					val = data[header->count-1];
					return true;// success
				} else {
					return false; // failure
				}
			}

			//get randomly: forbidden, please use the iterator!
			//
			data_iterator begin(){
				return &data[0];
			}
			data_iterator end(){
				return &data[header->count];
			}

			//Be careful to call!
			inline byte* address(){
				return (byte*)header;
			}
			inline Header_t* getHeader(){
				return header;
			}
		private:
			Header_t* header;
			Data_t* data;
			size_t size_;
		};

		template <typename Header_t, typename Data_t>
		class parFormatedBuffer
			:public FormatedBuffer<Header_t, Data_t> {
		private:
			typedef GRE::Util::TicketLock lock_t;
			typedef FormatedBuffer<Header_t, Data_t> Base;
			lock_t lock;
			bool isFull;
		public:
			parFormatedBuffer():lock(),isFull(false){}
			~parFormatedBuffer(){}
			inline bool push_back(const Data_t& val){
				bool ret;
				lock.spinLock();
				ret = Base::push_back(val);
				lock.unLock();
				if(!ret) isFull = true;
				return ret;
			}
			inline bool pop_back(Data_t& val){
				bool ret;
				lock.spinLock();
				ret = Base::pop_back(val);
				lock.unLock();
				return ret;
			}
			inline void pop_back(){
				lock.spinLock();
				Base::pop_back();
				lock.unLock();
			}
			inline bool back(Data_t& val){
				bool ret;
				lock.spinLock();
				ret = Base::back(val);
				lock.unLock();
				return ret;
			}
			//if return NULL, then it has been renewed by another thread.
			inline byte* renew(byte* addr, const size_t newSize){
				byte* ret = NULL;
				lock.spinLock();
				if(isFull == true){
					ret = Base::address();
					Base::format(addr, newSize);
					isFull = false;
				}
				lock.unLock();
				return ret;
			}
		};
	}//end Util
}

#endif
