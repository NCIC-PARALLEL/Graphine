#include <iostream>
#include <vector>
#include <string>

#include "GRE.hpp"
#include "Kernels/BC.hpp"

typedef GRE::Util::BufferPool BufferPool_t;
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;

typedef Forward Forward_t;
typedef Backward Backward_t;
typedef GRE::SynchronousEngine<Forward_t> ForwardEngine_t;
typedef GRE::SynchronousEngine<Backward_t> BackwardEngine_t;

int main(int argc, char** argv)
{	
	GRE::COMM::Underlying comm;
	comm.init(argc, argv);

	BufferPool_t pool;
	pool.config(2048, 8192*256, true);
	//pool.config(512, 512*256, true);

	std::string prefix(argv[1]);
	std::string format;
	std::string method;
	if(argc > 2){
		std::string tmp(argv[2]); 
		if(tmp=="bin")
			format = "bin";
		else
			format = "ascii";
	} else//in default
			format = "ascii";
	if(argc > 3){
		std::string tmp(argv[3]);
		if(tmp == "partition")
			method = "partition";
		else
			method = "edgelists";
	} else//in default 
		method = "edgelists";
	///////////////////////////
	DistributedGraph_t distributedGraph;
	comm.barrier();
	distributedGraph.ingress<InOut>(comm, pool, prefix, format, method);
	comm.barrier();
	//
	std::cout<<comm.rank()<<"Finishing loading...."<<std::endl;
	SourceVertexSet_t srcSet;
	std::vector<vid_t> src;
	const int n = GRE::Util::loadFromFile(prefix+".roots", src, 0);
	if(n==0){
		std::cerr<<"Err in reading src vertex from file."<<std::endl;
		exit(0);
	}
	srcSet.addVertex(src[1], distributedGraph);
	///////////////////////////

	///////////////////////////
	//phase-1
	vtxData* vtxState=NULL;//allocated in bfs
	Forward_t bfs(distributedGraph, vtxState);
	ForwardEngine_t bfsengine(distributedGraph, comm, pool, bfs);
	bfsengine.initialize();

	comm.initProfiling();
	comm.startTimer();
	bfsengine.run<Out>(srcSet, GRE::Par_ThreadPool);
	//bfsengine.run<Out>(srcSet, GRE::Par_ThreadGroup);
	double time1 = comm.stopTimer();
	//bfsengine.output_stat();
	comm.outputProfilingData();
	if(comm.rank()==0){
		std::cout<<"| Iteration count: "<<bfsengine.getSuperstepCount()<<std::endl;
		std::cout<<"| Use time: "<<time1<<"s"<<std::endl;
		std::cout<<comm.rank()<<"| init vertex time: "<<bfsengine.getInitVertexTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"| sync scatter time: "<<bfsengine.getSyncScatterTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"| main compute time: "<<bfsengine.getLocalComputeTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"| sync combiner time: "<<bfsengine.getSyncCombinerTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"| locla apply time: "<<bfsengine.getLocalApplyTime()<<"s"<<std::endl;
		//std::cout<<comm.rank()<<"|ThreadGroup mode | vertex access conflict rate: "<<engine.getVertexAccessConflictRate()<<std::endl;
	}
	int nLevels=bfs.getCurrLevel();
	bfsengine.finalize();
	///////////////////////////

	///////////////////////////
	//phase-2
	double* dependencies=NULL;
	//assert(vtxData);
	Backward_t bp(distributedGraph, vtxState, dependencies);
	BackwardEngine_t bpengine(distributedGraph, comm, pool, bp);
	bpengine.initialize();

	comm.startTimer();
	for(int i=nLevels-2; i>=0; i--){
		//std::cout<<"Backward Prop....."<<i<<std::endl;
		srcSet.setALL();
		bp.setCurrLevel(i);
		//bpengine.run<In>(srcSet, GRE::Par_ThreadGroup);
		bpengine.run<In>(srcSet, GRE::Par_ThreadPool);
	}
	double time2 = comm.stopTimer();
	if(comm.rank()==0){
		std::cout<<"ThreadPool mode | Use time: "<<time2<<"s"<<std::endl;
	}
	bpengine.finalize();
	///////////////////////////
	delete []vtxState;
	delete []dependencies;
	comm.finalize();
	pool.finalize();
}
