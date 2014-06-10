#ifndef GRE_Util_Locks_HPP
#define GRE_Util_Locks_HPP
#include "Types.hpp"
#include "Util/Atomic.hpp"
namespace GRE{
	namespace Util{
	class TicketLock
	{
	public:
		static const int size = 4;
		TicketLock()
		{
			lock[0].u = 0;
		}
		~TicketLock(){}

		inline void init(){
			lock[0].u = 0;
		}
		inline void reset(){
			lock[0].u = 0;
		}

		inline bool tryLock()
		{
			//!!! Don't use it! It is not a real tryLock as you think.
			spinLock();
			return true;
		}
		inline void spinLock()
		{
			uint16_t me = uint16_fetch_add( &(lock[0].s.users), 1);
			while ( lock[0].s.ticket != me ) ; //cpu_relax();
		}
		inline void unLock()
		{
			sync();
			lock[0].s.ticket++;
		}
	private:
		union TicketLock_
		{
			uint32_t u;
			struct
			{
				uint16_t ticket;
				uint16_t users;
			} s;
		};
	private:
		TicketLock_ lock[1];
	};

	class xcasLock
	{
	public:
		static const int size = 4;
		xcasLock():lock(0){}
		xcasLock(bool isLocked)
		{
			if(isLocked)
				lock = 1;
			else
				lock = 0;
		}
		~xcasLock(){}

		inline bool tryLock()
		{
			if(lock == 0)
				return (int32_cas(&lock, 0, 1));
			return false;
		}
		inline void spinLock()
		{
			while(!int32_cas(&lock, 0, 1))
			{
				while(lock == 1);
			}
		}
		inline void unLock()
		{
			sync();
			lock = 0;
		}
		inline void init()
		{
			lock = 0;
		}
		inline void reset()
		{
			lock = 0;
		}
	private:
		volatile int32_t lock;
	};
	}//end Util
}
#endif
