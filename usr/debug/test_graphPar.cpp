#include <iostream>
#include <vector>
#include <string>

#include "Types.hpp"
#include "Graph/GraphIngress.hpp"

typedef GRE::Graph::EdgePartition EdgePartition_t;

int main(int argc, char** argv)
{	
	if(argc < 2){
		std::cerr<<"Usage: exe graphName nParts"<<std::endl;
		exit(-1);
	}
	std::string prefix(argv[1]);
	const int n = atoi(argv[2]);
	///////////////////////////
	EdgePartition_t par;
	par.divideEdges(prefix, n);
	///////////////////////////
}
