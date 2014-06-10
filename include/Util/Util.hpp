#ifndef GRE_Util_Util_HPP
#define GRE_Util_Util_HPP

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "Util/Atomic.hpp"
#include "Types.hpp"
namespace GRE{
	namespace Util{
	//can be parallelized.
	inline int count1s(const unsigned long word){
		int n = 0;
		unsigned long bits = word;
		while(bits){
			n++;
			bits = bits & (bits-1);
		}
		return n;
	}

	long count1s(unsigned long* bitmap, const long size)
	{
		long n = 0;
		#pragma omp parallel
		{
			long count = 0;
			#pragma omp for
			for(long i = 0; i<size; i++){
				count+=count1s(bitmap[i]);
			}
			AtomicFetchAdd(&n, count);
		}
		return n;
	}

	template <typename T>
	void BitMap2Set(unsigned long* bitmap, const long size, std::vector<T>& keys){
		const long nKeys = count1s(bitmap, size);
		keys.resize(nKeys);
		long n = 0;
		for(long i = 0; i < size; i++){
			T key = i*sizeof(unsigned long);
			unsigned long tmp = bitmap[i];
			while(tmp){
				if(tmp & 0x1UL)
					keys[n++] = key;
				tmp >>=1;
				key++;
			}
		}
		//assert(n == nKeys);
	}

	template <typename T>
	inline T getMax(std::vector<T> vec){
		T max;
		typename std::vector<T>::iterator it = vec.begin();
		typename std::vector<T>::iterator end = vec.end();
		if(it!=end) max = *it;
		for(; it!=end; it++){
			if(*it > max)
				max = *it;
		}
		return max;
	}
	template <typename T>
	inline T getMin(std::vector<T> vec){
		T min;
		typename std::vector<T>::iterator it = vec.begin();
		typename std::vector<T>::iterator end = vec.end();
		if(it!=end) min = *it;
		for(; it!=end; it++){
			if(*it < min)
				min = *it;
		}
		return min;
	}


	template <typename T>
	T min(const T& t1, const T& t2){
		return (t1<t2? t1:t2);
	}
	template <typename T>
	T max(const T& t1, const T& t2){
		return (t1>t2? t1:t2);
	}

	template <typename T>
	std::string toString(T& val){
		std::stringstream ss;
		ss<<val;
		return ss.str();
	}

	struct Word{
		union {
			uint64_t word;
			uint32_t ints[2];
		};
		Word(const uint64_t inword):word(inword){}
	};
	template <typename T>
	int interpretWord(T* buf, const uint64_t bits, const int base){
		Word w(bits);
		int n=0;
		for(int i=0; w.word!=0; i++){
			uint32_t tmp=w.ints[i];
			int ind = base + 32*i;
			for(int j=0; tmp!=0; j++, tmp>>=1){
				if(tmp & 0x1){
					buf[n++]= ind + j;
				}
			}	
			w.ints[i]=0;
		}
		return n;
	}
	template <typename T>
	T align(T value, const int size = 64){
		if(value % size == 0)
			return value;
		else
			return (value+size-(value%size));
	}
	template <typename T>
	T uppAlign(T value){
		const T ret = value + (64 - value%64);
		return ret;
	}
	template <typename T>
	T downAlign(T value){
		const T ret = value - (value%64);
		return ret;
	}

	//
	template<typename T>
	size_t loadFromFile(const std::string prefix, std::vector<T>& array, const size_t n = 0){
		std::ifstream infile(prefix.c_str());
		if(!infile){
			std::cerr<<"Err: fail to read data from"<<prefix<<std::endl;
			return 0;
		}
		size_t count = 0;
		while(infile.good() && !infile.eof()){
			T elem;
			infile>>elem;
			array.push_back(elem);
			count++;
			if(n!=0  && count==n)
				break;
		}
		return count;
	}
	}//end Util
}
#endif
