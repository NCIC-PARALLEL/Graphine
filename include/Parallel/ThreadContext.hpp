#ifndef _GRE_Engine_Context_HPP
#define _GRE_Engine_Context_HPP
#include "Types.hpp"
#include "Parallel/ThreadGroup.hpp"

namespace GRE {
	template <typename T, typename VertexSet_t, typename Engine_t>
	class ThreadContext{
	public:
		typedef typename T::vLock_t vLock_t;
		typedef typename vLock_t::lock_t Lock_t;
	public:
		ThreadContext(T& _t, VertexSet_t& _nextSet, Engine_t* _engine)
			:nGroups(0),threadID(0),groupID(0),vlock(_t.refvLock()),nextSet(_nextSet),engine(_engine),n(0), c(0) {}

		~ThreadContext(){}
	public:
		//1
		int threadID;
		int nGroups;
		int groupID;
		//2
		VertexSet_t& nextSet;
		vLock_t& vlock;
		Engine_t* engine;
		//3:debug-only
		long n;
		long c;
	};

	template <typename VertexSet_t, typename Engine_t>
	class TG_ThreadContext{
	public:
		typedef GRE::ThreadGroup ThreadGroup_t;
		typedef typename ThreadGroup_t::Thread_t Thread_t;
		typedef typename ThreadGroup_t::vLock_t vLock_t;
		typedef typename vLock_t::lock_t Lock_t;
		typedef typename Engine_t::SocketFormatedBuffer_t SocketFormatedBuffer_t;
	public:
		TG_ThreadContext(Thread_t& _t, VertexSet_t& _nextSet, Engine_t* _engine)
			:group(_t.refThreadGroup()),threadID(_t.id),
			groupID(_t.refThreadGroup().getGroupID()),
			nGroups(_t.refThreadGroupShop().getNumGroups()),
			g_vlock(_t.refThreadGroup().refvLock()),
			vlock(_t.refThreadGroupShop().refvLock()),
			nextSet(_nextSet),engine(_engine),
			n(0), c(0){
			socketBuffers.resize(nGroups);
		}
		~TG_ThreadContext(){}
	public:
		ThreadGroup_t& group;
		//1
		const int threadID;
		const int groupID;
		const int nGroups;
		//2
		vLock_t& g_vlock;
		vLock_t& vlock;
		VertexSet_t& nextSet;
		Engine_t* engine;
		//3
		std::vector<SocketFormatedBuffer_t> socketBuffers;
		//
		long n;
		long c;
	};
}
#endif
