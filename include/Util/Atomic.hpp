#ifndef GRE_ATOMIC_HPP
#define GRE_ATOMIC_HPP
//#if defined(GCC_LINUX)
#if 1
extern "C"
{
#include <stdlib.h>
//wraped atomic operations, attention: no variable type check!!!
//auto-reuse
/*
short atomic_fetch_add(short* u,short v) 
{
	return __sync_fetch_and_add((u),(v));
}
short atomic_add_fetch(short* u,short v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(short* u,short v,short w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}

unsigned short atomic_fetch_add(unsigned short* u,unsigned short v) 
{
	return __sync_fetch_and_add((u),(v));
}
unsigned short atomic_add_fetch(unsigned short* u,unsigned short v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(unsigned short* u,unsigned short v,unsigned short w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}

int atomic_fetch_add(int* u,int v) 
{
	return __sync_fetch_and_add((u),(v));
}
int atomic_add_fetch(int* u,int v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(int* u,int v,int w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}

unsigned int atomic_fetch_add(unsigned int* u,unsigned int v) 
{
	return __sync_fetch_and_add((u),(v));
}
unsigned int atomic_add_fetch(unsigned int* u,unsigned int v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(unsigned int* u,unsigned int v,unsigned int w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}

long atomic_fetch_add(long* u,long v) 
{
	return __sync_fetch_and_add((u),(v));
}
long atomic_add_fetch(long* u,long v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(long* u,long v,long w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}

unsigned long atomic_fetch_add(unsigned long* u,unsigned long v) 
{
	return __sync_fetch_and_add((u),(v));
}
unsigned long atomic_add_fetch(unsigned long* u,unsigned long v) 
{
	return __sync_add_and_fetch((u),(v));
}
unsigned long atomic_fetch_or(unsigned long* u,unsigned long v) 
{
	return __sync_add_and_fetch((u),(v));
}
unsigned long atomic_or_fetch(unsigned long* u,unsigned long v) 
{
	return __sync_add_and_fetch((u),(v));
}
bool atomic_cas(unsigned long* u,unsigned long v,unsigned long w) 
{
	return __sync_bool_compare_and_swap((u),(v),(w));
}
*/
//explicit typenaming
//
#define AtomicFetchAdd(u,v) __sync_fetch_and_add((u),(v))
#define AtomicAddFetch(u,v) __sync_add_and_fetch((u),(v))
#define AtomicFetchOr(u,v) __sync_fetch_and_or((u),(v))
#define AtomicOrFetch(u,v) __sync_or_and_fetch((u),(v))
#define AtomicFetchAnd(u,v) __sync_fetch_and_and((u),(v))
#define AtomicAndFetch(u,v) __sync_and_and_fetch((u),(v))
#define AtomicCAS(u,v,w)	__sync_val_compare_and_swap((u),(v),(w))
#define sync() __sync_synchronize()

#define int64_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define int64_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define int64_cas(u,v,w) __sync_val_compare_and_swap((u),(v),(w))
#define uint64_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define uint64_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define uint64_cas(u,v,w) __sync_val_compare_and_swap((u),(v),(w))
#define int32_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define int32_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define int32_cas(u,v,w) __sync_bool_compare_and_swap((u),(v),(w))
#define uint32_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define uint32_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define uint32_cas(u,v,w) __sync_bool_compare_and_swap((u),(v),(w))
#define int16_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define int16_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define int16_cas(u,v,w) __sync_bool_compare_and_swap((u),(v),(w))
#define uint16_fetch_add(u,v) __sync_fetch_and_add((u),(v))
#define uint16_add_fetch(u,v) __sync_add_and_fetch((u),(v))
#define uint16_cas(u,v,w) __sync_bool_compare_and_swap((u),(v),(w))

}
#endif
#endif
