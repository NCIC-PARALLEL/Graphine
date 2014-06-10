#ifndef GRE_fineLock_HPP
#define GRE_fineLock_HPP
#include "Types.hpp"
#include "Atomic.hpp"
#include <cassert>
#include <cstdlib>

namespace GRE{
	namespace Util{
	class fLock_ID{
	public:
		typedef int32_t lock_t;	
	private:
		lock_t* locks;
	private:
		//to forbid implicit construction.
		fLock_ID(){}
	public:
		fLock_ID(const int n){
			assert(n>0);
			locks = (lock_t*)calloc(n, sizeof(lock_t));
			assert(locks);
		}
		~fLock_ID(){
			free(locks);
		}
		inline bool tryLock(const size_t vid){
			const size_t lid = vid;
			if(locks[lid])//locked
				return false;
			else
				return int32_cas(&locks[lid], 0, 1);
		}
		inline void spinLock(const size_t vid){
			const size_t lid = vid;
			while(locks[lid] || !int32_cas(&locks[lid], 0, 1));
		}
		inline void unLock(const size_t vid){
			const size_t lid = vid;
			sync();
			locks[lid] = 0;
		}
		inline lock_t& getLockRef(const size_t vid){
			const size_t lid = vid;
			return locks[lid];
		}
	};
	}//end Util
}
#endif
