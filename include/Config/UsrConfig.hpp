#ifndef GRE_USR_CONFIG_HPP
#define GRE_USR_CONFIG_HPP
#include "SysConfig.hpp"
#include "Types.hpp"
namespace GRE{
	//Buf Size
	const static size_t smallBufSize = pageSize*4;//4KB * 4 = 16KB
	const static size_t middleBufSize = pageSize*1024;//4KB * 1K = 4MB
	const static size_t largeBufSize = pageSize*1024*128;//4KB * 128K = 512MB
	const static size_t hugeBufSize = pageSize*1024*512;//4KB * 128K = 2GB

	const static size_t defaultBufPoolSize = 256;
	const static size_t defaultBufPoolSizeX1 = 1024;
	const static size_t defaultBufPoolSizeX2 = 4096;
	const static size_t defaultBufPoolSizeX3 = 4096*4;

	const static int MAX_NUM_MACHINES = 64;
	const static int MAX_NUM_THREADS = 32;
}
#endif
