#ifndef GRE_TIMER_HPP
#define GRE_TIMER_HPP
extern "C"
{
//This code fragmant was taken from graph500 benchmark, see more on http://www.graph500.org
#include "stdlib.h"
#include "unistd.h"
#include "stdlib.h"
#include "time.h"
#if defined(CLOCK_MONOTONIC)
#define TICTOC_CLOCK CLOCK_MONOTONIC
#define TICTOC_CLOCK_NAME "CLOCK_MONOTONIC"
#elif defined(CLOCK_REALTIME)
#define TICTOC_CLOCK CLOCK_REALTIME
#define TICTOC_CLOCK_NAME "CLOCK_REALTIME"
#else
#error "Failed to find a timing clock."
#endif
}
namespace GRE
{
class Timer
{
	public:
	Timer(){}
	~Timer(){}
	private:
	struct timespec tic_ts;
	public:
	void start()
	{
		clock_gettime(TICTOC_CLOCK, &tic_ts);	
	}
	double stop()
	{
		double interval;
		struct timespec ts;
		clock_gettime(TICTOC_CLOCK, &ts);
		interval = (ts.tv_nsec - (double)tic_ts.tv_nsec) * 1.0e-9;
		interval += (ts.tv_sec - (double)tic_ts.tv_sec);
		return interval;
	}
};
static unsigned long cpuHz=0;
class TSCTimer
{
public:
	TSCTimer():elapsedCycles(0),lastDida(0){}
	~TSCTimer(){}
	inline void init(){
		if(cpuHz == 0)
			cpuHz = calculateCPUhz();
		elapsedCycles = 0; 	
	}
	inline void start()
	{
		lastDida = read_tsc(); 	
	}
	inline void reset()
	{
		elapsedCycles = 0; 	
	}
	inline void stop()
	{
		elapsedCycles += (read_tsc() - lastDida); 	
	}
	inline double getTime()
	{
		if(cpuHz==0)
			cpuHz = calculateCPUhz();
		return (double)elapsedCycles/(double)cpuHz;
	}
	inline unsigned long getCycles(){
		return elapsedCycles;
	}
	inline unsigned long getCPUhz(){
		if(cpuHz==0)
			cpuHz = calculateCPUhz();
		return cpuHz;
	}

private:
	unsigned long elapsedCycles;
	unsigned long lastDida;
private:
	static inline unsigned long  calculateCPUhz()
	{
		unsigned long lastTemp=read_tsc();
		sleep(1);
		return (read_tsc()-lastTemp);
	}
	static inline unsigned long read_tsc()
 	{
 		unsigned int low,high;
 		unsigned long val;
 		__asm__ volatile("rdtsc" : "=a"(low),"=d"(high));
 		val=low|((unsigned long)(high)<<32);
 		return val;
 	}
};
//unsigned long TSCTimer::cpuHz = getCPUhz();
}
#endif
