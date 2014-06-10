#include <iostream>
#include <vector>
#include <string>

#include "GRE.hpp"
#include "Kernels/BFS.hpp"

typedef GRE::Util::BufferPool BufferPool_t;
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;
typedef BFS BFS_t;
typedef GRE::SynchronousEngine<BFS_t> Engine_t;

int main(int argc, char** argv)
{	
	GRE::COMM::Underlying comm;
	comm.init(argc, argv);

	BufferPool_t pool;
	pool.config(4096, 8192*256, true);
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
	distributedGraph.ingress<InOut>(comm, pool, prefix, format, method);
	comm.barrier();

	//distributedGraph.save(prefix);
	//Load source vertex set
	SourceVertexSet_t srcSet;
	//secSet.loadFromFile(prefix);
	//srcSet.setFULL();
	std::vector<vid_t> src;
	const int n = GRE::Util::loadFromFile(prefix+".roots", src, 0);
	if(n==0){
		std::cerr<<"Err in reading src vertex from file."<<std::endl;
		exit(0);
	}
	srcSet.addVertex(src[1], distributedGraph);

	BFS bfs(distributedGraph);
	Engine_t engine(distributedGraph, comm, pool, bfs);
	engine.initialize();
	//engine.configMaxSupersteps(1);
	//comm.barrier();

	//run with tp
	comm.initProfiling();
	comm.startTimer();
	engine.run<In>(srcSet, GRE::Par_ThreadPool);
	double time1 = comm.stopTimer();
	comm.outputProfilingData();
	if(comm.rank()==0){
		std::cout<<"ThreadPool mode | Iteration count: "<<engine.getSuperstepCount()<<"s"<<std::endl;
		std::cout<<"ThreadPool mode | Use time: "<<time1<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | init vertex time: "<<engine.getInitVertexTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync scatter time: "<<engine.getSyncScatterTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | main compute time: "<<engine.getLocalComputeTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync combiner time: "<<engine.getSyncCombinerTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | apply time: "<<engine.getLocalApplyTime()<<"s"<<std::endl;
	}
	///////////////////////////
	engine.finalize();
	comm.finalize();
	pool.finalize();
}
