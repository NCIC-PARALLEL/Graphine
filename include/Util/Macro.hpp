#ifndef _GRE_Util_Macro_HPP
#define _GRE_Util_Macro_HPP

//1:Compiler implications for branch predication.
//However, in our experiments, it seems that this trick makes no sense.
#undef likely
#undef unlikely

#if defined(__GNUC__) && __GNUC__ >= 4
#define likely(x)   (__builtin_expect((x), 1))
#define unlikely(x) (__builtin_expect((x), 0))
#else
#define likely(x)   (x)
#define unlikely(x) (x)
#endif

//2:Container_of
#define container_of(ptr, type, member) (type *)((char *)(ptr) - offsetof(type, member))

//3:

#endif
