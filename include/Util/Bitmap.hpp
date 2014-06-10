#ifndef GRE_Util_Bitmap_HPP
#define GRE_Util_Bitmap_HPP

#include "Types.hpp"
#include "Util/Atomic.hpp"
#include "Util/Prime.hpp"
#include <cstring>
#include <cstdlib>
namespace GRE
{
	namespace Util{
	class Bitmap
	{
	//one word is an unsigned long, i.e. one word == 64 bits.
	public:
		Bitmap():size_(0),bits(NULL){}
		explicit Bitmap(int64_t nBits)
		{
			if(nBits > 0){
				size_ = nBits;
				int64_t nWords = size_/64 + 1;
				bits = (uint64_t*)calloc(nWords, sizeof(uint64_t));
			}
		}
		~Bitmap()
		{
			if(bits) free(bits);
		}
		void destroy(){
			if(bits) {
				free(bits);
				bits = NULL;
			}
		}
		#ifdef Debug
		bool testValid(){
			assert(sizeof(uint64_t)==8);
		}
		#endif
		bool resize(const int64_t nBits)
		{
			if(nBits > size_){
				size_ = nBits;
				if(bits) free(bits);
				int64_t nWords = size_/64 + 1;
				bits = (uint64_t*)calloc(nWords, sizeof(uint64_t));
			}
			return (bits!=NULL);
		}
		inline void set(int64_t idx)
		{
			AtomicFetchOr(&bits[idx/64],0x1UL<<(idx%64));
		}
		inline bool test(int64_t idx)
		{
			//if(bits[idx/64] & (0x1UL<<(idx%64))) 
			//	return true;;
			//return false;;
			return (bits[idx/64] & (0x1UL<<(idx%64)));
		}
		inline bool testAndSet(int64_t idx)
		{
			if(bits[idx/64]&(1UL<<(idx%64)))
				return true;
			return ((AtomicFetchOr(&bits[idx/64],(1UL<<(idx%64)))) & (1UL<<(idx%64)));
		}
		inline void unset(int64_t idx)
		{
			const uint64_t val = ~(0x1L<<(idx%64));
			AtomicFetchAnd(&bits[idx/64], val);
		}
		inline int64_t size()
		{
			return size_;
		}
		void clear()
		{
			memset(bits, 0, size_/8+8);
		}
		void reset()
		{
			clear();
		}
	private:
		int64_t size_;
		uint64_t* bits;
	public:
		typedef uint64_t* wordIterator;
		wordIterator begin(){
			return &bits[0];
		}
		wordIterator end(){
			return &bits[size_/64+1];
		}
		wordIterator getWordIterator(int64_t id){
			return &bits[id/64];
		}

		inline void setWord(int64_t windex, const uint64_t newWord){
			bits[windex] = newWord;
		}
		inline uint64_t getWord(int64_t windex){
			return bits[windex];
		}
		inline int64_t nWords(){
			return (size_/64+1);
		}
	public://serial interface
		inline void set_serial(const int64_t idx)
		{
			bits[idx/64] |= (0x1UL<<(idx%64));
		}
		inline void unset_serial(const int64_t idx)
		{
			const uint64_t val = ~(0x1L<<(idx%64));
			bits[idx/64] &= val;
		}
		template <typename ostream_t>
		void out(ostream_t& os){
			os.write(reinterpret_cast<char*>(&size_), sizeof(size_));
			os.write(reinterpret_cast<char*>(bits), size_/8+8);
		}
		template <typename istream_t>
		void in(istream_t& is){
			int64_t n;
			is.read(reinterpret_cast<char*>(&n), sizeof(int64_t));
			resize(n);
			is.read(reinterpret_cast<char*>(bits), size_/8+8);
		}
	};

	template <size_t size_=64>
	class FixedBitmap
	{
	//one word is an unsigned long, i.e. one word == 64 bits.
	public:
		FixedBitmap(){
			//memset(bits, 0, size_/8);
			//PLEASE: explitly call this->clear
		}
		~FixedBitmap(){}

		#ifdef Debug
		bool testValid(){
			assert(size_%64==0);
			assert(sizeof(uint64_t)==8);
		}
		#endif

		inline void set(int64_t idx)
		{
			AtomicFetchOr(&bits[idx/64],0x1UL<<(idx%64));
		}

		inline bool test(int64_t idx)
		{
			return (bits[idx/64] & (0x1UL<<(idx%64)));
		}

		inline bool testAndSet(int64_t idx)
		{
			if(bits[idx/64]&(1UL<<(idx%64)))
				return true;
			return ((AtomicFetchOr(&bits[idx/64],(1UL<<(idx%64)))) & (1UL<<(idx%64)));
		}

		inline void unset(int64_t idx)
		{
			const uint64_t val = ~(0x1L<<(idx%64));
			AtomicFetchAnd(&bits[idx/64], val);
		}

		inline int64_t size()
		{
			return size_;
		}

		void clear()
		{
			memset(bits, 0, size_/8);
		}

		template<typename stream_t>
		size_t dump(stream_t& os){
			os.write(bits, size_/8);
		}

	private:
		uint64_t bits[size_/64];

	public://Word-level operarions
		typedef uint64_t* wordIterator;
		wordIterator begin(){
			return &bits[0];
		}

		wordIterator end(){
			return &bits[size_/64];
		}

		inline void setWord(int64_t windex, const uint64_t newWord){
			bits[windex] = newWord;
		}

		inline uint64_t getWord(int64_t windex){
			return bits[windex];
		}

		inline int64_t nWords(){
			return (size_/64);
		}
		
		template <typename T>
		int exportIndexValues(T* buf){
			int n=0;
			for(int i=0; i<size_/64; i++){
				uint64_t tmp = bits[i];
				int index = 64*i;
				while(tmp){
					if(tmp & 0x1UL){
						buf[n++]=index;
					}
					index++;
					tmp >>=1;
				}
			}
			return n;
		}
	public://serial interface
		inline void set_serial(const int64_t idx){
			bits[idx/64] |= (0x1UL<<(idx%64));
		}

		inline void unset_serial(const int64_t idx){
			const uint64_t val = ~(0x1L<<(idx%64));
			bits[idx/64] &= val;
		}
	};
	}//end Namepace-Util
}
#endif
