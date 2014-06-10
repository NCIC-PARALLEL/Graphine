#ifndef GRE_Util_HashIndex_HPP
#define GRE_Util_HashIndex_HPP

#include <omp.h>
#include "Util/Prime.hpp"
#include <assert.h>
namespace GRE{
	namespace Util{
		//This is actually a Cuckoo Map.
		template <typename key_t, typename value_t, typename offset_t=int>
		class HashIndex{
		public:
			typedef std::pair<key_t, value_t> entry_t;
		public:
			HashIndex(const key_t invalidV):illegalValue(invalidV),hashSize(0), offsets(),entries(){}
			~HashIndex(){}
			//
			value_t value(const key_t key){
				const long idx = hash(key);
				for(long i = offsets[idx]; i<offsets[idx+1]; i++){
					if(entries[i].first == key) 
						return entries[i].second;
					else {
						if(entries[i].first > key)
							break;
					}
				}
				return illegalValue;
			}
			//
			void buildFromArray(std::vector<key_t>& keys)
			{
				buildFromArray(keys, 0, keys.size());
			}
			void buildFromArray(std::vector<key_t>& keys, const long beginIdx)
			{
				const long nKeys = keys.size()-beginIdx;
				buildFromArray(keys, beginIdx, nKeys);
			}
			void buildFromArray(std::vector<key_t>& keys, const long beginIdx, const long _nKeys)
			{
				assert(beginIdx >= 0 && _nKeys >0);
				assert(beginIdx + _nKeys <= keys.size());
				//
				const long nKeys = _nKeys;
				hashSize = getPrime(nKeys);
				//
				//
				offsets.resize(hashSize+1,0);
				entries.resize(nKeys);
				//
				//
				#pragma omp parallel for
				for(long i = 0; i<nKeys; i++){
					const offset_t off = hash(keys[beginIdx+i]);
					AtomicFetchAdd(&offsets[off],1);
				}
				offset_t off = offsets[0];
				offsets[0] = 0;
				for(long i=1; i<=hashSize; i++){
					const offset_t off1 = off;
					off = offsets[i];
					offsets[i] = off1+offsets[i-1]; 
				}
				//
				//
				#pragma omp parallel
				{
					const int nThreads = omp_get_num_threads();
					const int tid = omp_get_thread_num();
					for(long i = 0; i<nKeys; i++){
						const long off = hash(keys[beginIdx+i]);
						if(off % nThreads == tid){
							entries[offsets[off]].first = keys[beginIdx+i];
							entries[offsets[off]].second = beginIdx+i;
							offsets[off]++;
						}
					}
				}	
				//
				//
				for(long i=hashSize; i>0; i--){
					offsets[i] = offsets[i-1]; 
				}
				offsets[0]=0;
			}

			void clear(){
				hashSize=0;
				offsets.clear();
				entries.clear();
			}
			void reset(){
				clear();
			}
		private:
			const key_t illegalValue;
			long hashSize;
			std::vector<offset_t> offsets;
			std::vector<entry_t> entries;
		private:
			inline long hash(const key_t key)
			{
				return key%hashSize;
			}
		public:
			//Binary IO Interface
			void load(){
				//load
			}
			void save(){
				//save
			}
		};
	}//end Util
}
#endif
