#ifndef GRE_SYSCONFIG_HPP
#define GRE_SYSCONFIG_HPP
namespace GRE{
	const static size_t pageSize = 4096; 
	const static int cacheLineSize = 64;//bytes
	//
	const static int Sys_nSockets = 2;
	const static int Sys_nCoresPerSocket = 6;
	const static int Sys_nThreadsPerCore = 1;
}
#endif
