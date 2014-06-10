#ifndef GRE_COMM_Format_HPP
#define GRE_COMM_Format_HPP
#include "Types.hpp"
//#include <stdint>
namespace GRE{
	namespace COMM{
	class Protocol{
	public:
		typedef byte opcode_t;
		typedef byte flag_t;

		static const opcode_t op_NOP = 0;
		static const opcode_t op_END = 1;
		static const opcode_t op_invalid = 255;//e.g. -1
		//
		//definition of flag
		static const flag_t flag_end = 1;
		static const flag_t flag_reponse = 2;
		static const flag_t flag_noresponse = 4;
		static const flag_t flag_4 = 8;
		static const flag_t flag_5 = 16;
		static const flag_t flag_6 = 32;
		static const flag_t flag_7 = 64;

		union Header{
				uint64_t u;
				struct{
					opcode_t op;
					flag_t flag;
					uint16_t extra;//e.g. dest
					int32_t count;
				}; 
			Header():u(0){}
			~Header(){}
		public:
			//definition of Op: [0, 255]
			//Public Interface
			inline void clear(){
				u = 0;
			}
			inline void setOp(opcode_t opcode){
				op = opcode;
			}
			inline opcode_t getOp(){
				return op;
			}
			inline void setFlag(flag_t iflag, const flag_t imask = -1){
				flag &= imask;
				flag |= iflag;
			}
			inline bool testFlag(flag_t iflag){
				return ((flag & iflag)!=0);
			}
			inline flag_t getFlag(){
				return flag;
			}
			inline void setExtra(uint16_t iextra){
				extra = iextra;
			}
			inline uint16_t getExtra(){
				return extra;
			}
			inline void setCount(uint32_t icount){
				count = icount;
			}
			inline uint32_t getCount(){
				return count;
			}
			//begin debug: should be removed later!!
			inline bool isEnd(){
				return testFlag(flag_end);//255
			}
			inline void setEnd(){
				setFlag(flag_end);
			}
			//end debug
		};
	};
	}//end COMM
}
#endif
