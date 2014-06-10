#ifndef GRE_vLock_HPP
#define GRE_vLock_HPP
#include "Types.hpp"
#include "Atomic.hpp"
#include "Util/Locks.hpp"

#define lockSize 8192L
#define lockSizePrime 8191UL
#define lockMask 8191UL

//#define USE_Directly
namespace GRE{
	namespace Util{
#ifdef USE_Directly
	class vLock_ID{
	public:
		typedef int32_t lock_t;	
	private:
		lock_t locks[lockSize];
	public:
		vLock_ID(){}
		~vLock_ID(){};
		inline bool tryLock(const size_t vid){
			const size_t lid = vid & lockMask;
			if(locks[lid])//locked
				return false;
			else
				return int32_cas(&locks[lid], 0, 1);
		}
		inline void spinLock(const size_t vid){
			const size_t lid = vid & lockMask;
			while(locks[lid] || !int32_cas(&locks[lid], 0, 1));
		}
		inline void unLock(const size_t vid){
			const size_t lid = vid & lockMask;
			sync();
			locks[lid] = 0;
		}
		inline lock_t& getLockRef(const size_t vid){
			const size_t lid = vid & lockMask;
			return locks[lid];
		}
	};

	class vLock_AD{
	public:
		typedef int32_t lock_t;	
	private:
		lock_t locks[lockSize];
	public:
		vLock_AD(){}
		~vLock_AD(){};
		inline bool tryLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
				if(locks[lid])//locked
				return false;
			else
				return int32_cas(&locks[lid], 0, 1);
		}
		inline void spinLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			while(locks[lid] || !int32_cas(&locks[lid], 0, 1));
		}
		inline void unLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			sync();
			locks[lid] = 0;
		}
		inline lock_t& getLockRef(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			return locks[lid];
		}
	};
#else
	class vLock_ID{
	public:
		typedef GRE::Util::xcasLock lock_t;	
	private:
		lock_t locks[lockSize];
	public:
		vLock_ID(){}
		~vLock_ID(){};
		inline bool tryLock(const size_t vid){
			const size_t lid = vid & lockMask;
			return locks[lid].tryLock();
		}
		inline void spinLock(const size_t vid){
			const size_t lid = vid & lockMask;
			locks[lid].spinLock();
		}
		inline void unLock(const size_t vid){
			const size_t lid = vid & lockMask;
			locks[lid].unLock();
		}
		inline lock_t& getLockRef(const size_t vid){
			const size_t lid = vid & lockMask;
			return locks[lid];
		}
	};

	class vLock_AD{
	public:
		typedef GRE::Util::xcasLock lock_t;	
	private:
		lock_t locks[lockSize];
	public:
		vLock_AD(){}
		~vLock_AD(){};
		inline bool tryLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			return locks[lid].tryLock();
		}
		inline void spinLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			locks[lid].spinLock();
		}
		inline void unLock(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			locks[lid].unLock();
		}
		inline lock_t& getLockRef(void* add){
			const size_t lid = reinterpret_cast<unsigned long>(add) % lockMask;
			return locks[lid];
		}
	};
#endif
	}//end Util
}
#endif
