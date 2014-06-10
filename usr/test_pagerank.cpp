#include <iostream>
#include <vector>
#include <string>

#include "Types.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"
#include "Graph/GraphIngress.hpp"
#include "Graph/VertexSet.hpp"

#include "Comm/CommUnderlying.hpp"
#include "Comm/BufferingComm.hpp"
#include "Util/BufferPool.hpp"
#include "Engine/Engine.hpp"
//#include "Engine/EngineProfiling.hpp"
#include "Kernels/PageRank.hpp"

typedef GRE::Util::BufferPool BufferPool_t;
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;

//typedef PageRank_delta PageRank_t;
typedef PageRank PageRank_t;
typedef GRE::SynchronousEngine<PageRank_t> Engine_t;

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
	}
	///////////////////////////
	DistributedGraph_t distributedGraph;
	comm.barrier();
	distributedGraph.ingress(comm, pool, prefix, format, method);
	comm.barrier();

	//pool.finalize();
	//pool.config(4096, 8192*64, true);
	//distributedGraph.save(prefix);
	//Load source vertex set
	SourceVertexSet_t srcSet;
	//secSet.loadFromFile(prefix);
	srcSet.setALL();
	///////////////////////////
	PageRank_t pr(distributedGraph);
	Engine_t engine(distributedGraph, comm, pool, pr);
	engine.initialize();
	engine.configMaxSupersteps(10);
	//run with tg
	comm.initProfiling();
	comm.startTimer();
	engine.run(srcSet, GRE::Par_ThreadPool);
	//engine.run(srcSet,GRE::Par_ThreadGroup);
	double time1 = comm.stopTimer();
	//engine.output_stat();
	comm.outputProfilingData();
	if(comm.rank()==0){
		std::cout<<"ThreadPool mode | Iteration count: "<<engine.getSuperstepCount()<<"s"<<std::endl;
		std::cout<<"ThreadPool mode | Use time: "<<time1<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | init vertex time: "<<engine.getInitVertexTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync scatter time: "<<engine.getSyncScatterTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | main compute time: "<<engine.getLocalComputeTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync combiner time: "<<engine.getSyncCombinerTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | locla apply time: "<<engine.getLocalApplyTime()<<"s"<<std::endl;
		//std::cout<<comm.rank()<<"|ThreadGroup mode | vertex access conflict rate: "<<engine.getVertexAccessConflictRate()<<std::endl;
	}
	///////////////////////////
	engine.finalize();
	comm.finalize();
	pool.finalize();
}
