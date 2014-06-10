#ifndef GRE_Util_SparseIndex_HPP
#define GRE_Util_SparseIndex_HPP
#include "CodeSpec.hpp"
#include <omp.h>
namespace GRE{
	namespace Util{
		//
		class SparseIndex{
		public:
			typedef int offset_t;
			struct entry{
				unsigned long bitmap;
				offset_t offset;
			};
		private:
			vector<entry> L1Index;
			vector<entry> L2Index;
			vector<lvid_t> lvids;
		public:
			SparseIndex():L1Index(),L2Index(),lvids(){}
			~SparseIndex(){}
			lvid_t lvid(vid_t _vid){
				const unsigned long vid = getIndex(_vid);
				const unsigned long l2 = (vid >> 6) & 0x3fUL;
				const unsigned long l1 = vid >> 12;
				//locate value
				unsigned long bits = L1Index[l1].bitmap;
				if(bits & (0x1UL << l2)){
					bits = bits & ((0x1UL << 12) - 0x1UL);
					offset_t l2Offset = L1Index[l1].offset + count1s(bits);	
					unsigned long bits2 = L2Index[l2Offset].bitmap;
					const unsigned long l3 = vid & 0x3fUL;
					if(bits2 & (0x1UL << l3)){
						offset_t l3Offset = L2Index[l2Offset].offset + count1s(bits2);
						return lvids[l3Offset];
					}
				}
				return NotALvid;
			}
			void build(unsigned long* bitmap, long count){
				long n = count-1;
				while(!bitmap[n]) n--;
				vid_t max = n*sizeof(unsigned long);
				unsigned long tmp = bitmap[n];
				while(tmp>>=1) max++;

				lvids.reserve(max);
				L1Index.reserve(max>>12);
				//1:
				offset_t off = 0;	
				const long maxL1Index = max >> 12;
				for(long i=0; i<maxL1Index; i++){
					L1Index[i].offset = off;				
					L1Index[i].bitmap = 0;				
					for(long j=0; j<64; j++){
						if(bitmap[i*64+j]){
							L2Index.pushback(entry(bitmap[i*64+j], 0));				
							L1Index[i].bitmap |= (0x1UL << j);				
							off++;
						}
					}
				}
				//2:
				#pragma omp parallel for
				for(long i=0; i<L2Index.size(); i++){
					const int c = count1s(L2Index[i].bitmap);
					L2Index[i].offset = c;
				}
				//3:
				off = L2Index[0].offset;
				L2Index[0].offset = 0;
				for(long i=1; i<L2Index.size();i++)
				{
					const offset_t off1 = off;
					off = L2Index[i].offset;
					L2Index[i].offset = off1 + L2Index[i-1].offset;
				}
			}
		private:
			inline int count1s(unsigned long bits){
				int n = 0;
				while(bits){
					bits = bits & (bits-1);
					n++;
				}
				return n;
			}
		};
	}//end Util
}
#endif
