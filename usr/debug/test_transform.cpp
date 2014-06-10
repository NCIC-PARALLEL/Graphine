#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include "Util/Util.hpp"
int main(int argc, char** argv){
	std::string prefix(argv[1]);
	const int n = atoi(argv[2]);
	std::vector<size_t> longs;
	std::vector<int> ints;
	for(int i = 0; i < n; i++){
		const std::string fname1 = prefix + ".outdegree." + GRE::Util::toString(n) + "." + GRE::Util::toString(i);
		std::ifstream infile;
		infile.open(fname1.c_str(), std::ifstream::in | std::ifstream::binary);
		infile.seekg(0, std::ifstream::end);
		const long nV = (size_t)infile.tellg()/sizeof(size_t);
		std::cout << "trans.." << fname1 << "size.." << nV << std::endl << std::flush;
		longs.resize(nV);
		infile.read(reinterpret_cast<char*>(&longs[0]), nV*sizeof(size_t));
		infile.close();

		ints.resize(nV);
 		#pragma omp parallel for
		for(int j=0; j<nV; j++){
			ints[j] = longs[j];
		}
		std::ofstream outfile;
		outfile.open(fname1.c_str(), std::ofstream::out | std::ofstream::binary);
		outfile.write(reinterpret_cast<char*>(&ints[0]), nV*sizeof(int));
		outfile.close();
	}
}
