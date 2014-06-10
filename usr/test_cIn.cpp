#include <iostream>
#include <vector>
#include <string>

#include "GRE.hpp"
#include "Kernels/CountIn.hpp"

typedef GRE::Util::BufferPool BufferPool_t;
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;
typedef CountIn CountIn_t;
typedef GRE::SynchronousEngine<CountIn_t> Engine_t;

int main(int argc, char** argv)
{	
	GRE::COMM::Underlying comm;
	comm.init(argc, argv);

	BufferPool_t pool;
	pool.config(4096, 8192*128, true);

	std::string prefix(argv[1]);
	///////////////////////////
	DistributedGraph_t distributedGraph;
	comm.barrier();
	distributedGraph.ingress<InOut>(comm, pool, prefix);
	comm.barrier();
	//distributedGraph.save(prefix);
	//Load source vertex set
	SourceVertexSet_t srcSet;
	//secSet.loadFromFile(prefix);
	srcSet.setALL();
	///////////////////////////
	CountIn_t countIn(distributedGraph);
	Engine_t engine(distributedGraph, comm, pool, countIn);
	engine.initialize();
	engine.configMaxSupersteps(1);
	comm.startTimer();

	//comm.barrier();
	//engine.reset();

	comm.startTimer();
	engine.run<In>(srcSet, GRE::Par_ThreadPool);
	//engine.run(srcSet,GRE::Par_ThreadGroup);
	double time1 = comm.stopTimer();
	if(comm.rank()==0){
		std::cout<<"ThreadPool mode | Iteration count: "<<engine.getSuperstepCount()<<std::endl;
		std::cout<<"ThreadPool mode | Use time: "<<time1<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | init vertex time: "<<engine.getInitVertexTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync scatter time: "<<engine.getSyncScatterTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | main compute time: "<<engine.getLocalComputeTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | sync combiner time: "<<engine.getSyncCombinerTime()<<"s"<<std::endl;
		std::cout<<comm.rank()<<"|ThreadPool mode | locla apply time: "<<engine.getLocalApplyTime()<<"s"<<std::endl;
	}

	///////////////////////////
	engine.finalize();
	comm.finalize();
	pool.finalize();
}
