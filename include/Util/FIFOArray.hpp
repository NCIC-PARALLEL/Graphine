#ifndef GRE_Util_FIFOArray_HPP
#define GRE_Util_FIFOArray_HPP
#include <stdint.h>
namespace GRE{
	namespace Util{
		//Note: only support serial operations.
		template <typename T, int n=1024>
		class FIFOArray{
		public:
			FIFOArray():head(0),tail(0){}
			~FIFOArray(){}
			//interface
			inline bool get(T& val){
				if(head > tail){
					val = data[tail%n];
					tail++;
					return true;
				}
				return false;
			}
			//Note: Not safe, you have to call empty() to ensure the queue is not empty!
			inline T get(){
				const size_t tail1= tail;
				tail++;
				return data[tail1 % n];
			}
			inline bool put(T& val){
				if(head - tail < n){
					data[head % n] = val;
					head++;
					return true;
				}
				return false;
			}
			inline bool empty()const{
				return (head <= tail);
			}
			inline bool full()const{
				return (head - tail >= n);
			}
			void reset(){
				head = 0;
				tail = 0;
			}
			size_t count()const{
				return (head - tail); 
			} 	
			//TODO:debug only, delete later!!!
			inline size_t getTail(){
				return tail;
			}
		private:
			volatile size_t head;
			volatile size_t tail;
			T data[n];
		};
	}//end Util
}
#endif
