#include <string>
#include <iostream>
#include <fstream>
int main(int argc, char** argv){
	std::string prefix(argv[1]);
	std::string filename = prefix + ".edges";
	std::string filename1 = prefix + ".tsv";
	std::ifstream infile(filename.c_str());
	std::ofstream outfile(filename1.c_str());
	std::string meta;
	std::getline(infile, meta);
	std::cout<<meta<<std::endl;
	long cc = 0;
	while(infile.good() && !infile.eof()){
		cc++;
		std::string line;
		std::getline(infile, line);
		outfile<<line<<std::endl;
	}
	std::cout<<"ending..."<<cc<<std::endl;
	return 0;
}
