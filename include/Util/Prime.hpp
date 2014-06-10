#ifndef GRE_Util_Prime_HPP
#define GRE_Util_Prime_HPP
#include <cmath>
namespace GRE{
	namespace Util{
		template <typename T>
		T getPrime(T _value)
		{
			T candidate = _value;
			if(candidate <= 2) return 2;
			if(candidate == 3) return 3;
			if(candidate%2 == 0) candidate++; 
			while(true)
			{
				if(isPrime(candidate))
					return candidate;
				candidate+=2;
			}
		}
		template <typename T>
		static inline bool isPrime(T _cand)
		{
			if(_cand%2 == 0 || _cand%3 == 0) return false;
			T end = (T)sqrt((double)_cand);
			for(T i=5; i<=end; i+=2)
			{
				if((_cand % i) == 0)
					return false;
			}
			return true;
		}
	}//end Util
}
#endif
