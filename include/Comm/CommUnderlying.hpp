#ifndef GRE_COMM_Underlying1_HPP
#define GRE_COMM_Underlying1_HPP

//Note: Use the C-version of MPI routines!
//Don't distinguish data types, treat all as bytes(char).
#include <mpi.h>
#include "Types.hpp"
#include <vector>
namespace GRE{
	namespace COMM{
	class Underlying{
	public:
		typedef MPI_Op Op_t;
		typedef MPI_Datatype DataType_t;
	private:
		//Sata
		bool initialized;
		int _size;//number of processes
		int _rank;//my process id
		MPI_Comm comm;
		std::vector<MPI_Request> sendRequests;
		std::vector<MPI_Request> recvRequests;
		//Time measure
		double time;
		//Profiling
		#ifdef Comm_Profiling_Send
		long sendDataVolume;
		#endif
		#ifdef Comm_Profiling_Recv
		long recvDataVolume;
		#endif
	public:
		Underlying(const MPI_Comm _comm=MPI_COMM_WORLD):comm(_comm), initialized(false), time(0){}
		~Underlying(){}
		//Interface
		void init(int& argc,char**& argv )
		{
			if(!initialized){
				MPI_Init(&argc, &argv);
				MPI_Comm_size(comm, &_size);
				MPI_Comm_rank(comm, &_rank);
				sendRequests.resize(_size);
				recvRequests.resize(_size);
				#ifdef Comm_Profiling_Send
				sendDataVolume = 0;
				#endif
				#ifdef Comm_Profiling_Recv
				recvDataVolume = 0;
				#endif
				initialized = true;
			}
		}
		void finalize()
		{
			if(initialized){
				sendRequests.clear();
				recvRequests.clear();
				MPI_Finalize();
				initialized = false;
			}
		}
		//helpers
		inline int rank() const
		{
			return _rank;
		}
		inline int size() const
		{
			return _size;
		}
		//Blocking-Communication
		bool send(void* buf, int count, int dest, int tag = 0) //const
		{
			#ifdef Comm_Profiling_Send
			sendDataVolume += count;
			#endif
			return	(MPI_SUCCESS==MPI_Send(buf, count, MPI_BYTE, dest, tag, comm));
		}
		bool recv(void* buf, int count, const int src = MPI_ANY_SOURCE, const int tag = MPI_ANY_TAG) //const
		{
			#ifdef Comm_Profiling_Recv
			recvDataVolume += count;
			#endif
			return	(MPI_SUCCESS==MPI_Recv(buf, count, MPI_BYTE, src, tag, comm, MPI_STATUS_IGNORE));
		}
		//NonBlocking
		bool isend(void* buf, int count, int dest, int tag = 0)
		{
			#ifdef Comm_Profiling_Send
			sendDataVolume += count;
			#endif
			return (MPI_SUCCESS==MPI_Isend(buf, count, MPI_BYTE, dest, tag, comm, &sendRequests[dest]));	
		}
		bool irecv(void* buf, int count, const int src = MPI_ANY_SOURCE, const int tag = MPI_ANY_TAG)
		{
			#ifdef Comm_Profiling_Recv
			recvDataVolume += count;
			#endif
			//there is a trick here.
			if(src == MPI_ANY_SOURCE) 
				return (MPI_SUCCESS==MPI_Irecv(buf, count, MPI_BYTE, src, tag, comm, &recvRequests[_rank]));
			else
				return (MPI_SUCCESS==MPI_Irecv(buf, count, MPI_BYTE, src, tag, comm, &recvRequests[src]));
		}

		inline bool testRecv(int src){
			int flag;
			MPI_Test(&recvRequests[src], &flag, MPI_STATUS_IGNORE);
			if(flag)
				return true;
			else
				return false;
		}
		inline bool testRecv(){
			int flag;
			MPI_Test(&recvRequests[_rank], &flag, MPI_STATUS_IGNORE);
			if(flag)
				return true;
			else
				return false;
		}
		//
		inline bool testRecvExtraAnonymous(int* src, int* count, int* tag)
		{
			int flag;
			MPI_Status status;
			MPI_Test(&recvRequests[_rank], &flag, &status);
			if(flag){
				*src = status.MPI_SOURCE;
				*tag = status.MPI_TAG;
				MPI_Get_count(&status, MPI_BYTE, count);
				return true;
			}
			return false;
		}
		inline bool testRecvExtraTagAndCount(int src, int* count, int* tag)
		{
			int flag;
			MPI_Status status;
			MPI_Test(&recvRequests[src], &flag, &status);
			if(flag){
				*tag = status.MPI_TAG;
				MPI_Get_count(&status, MPI_BYTE, count);
				return true;
			}
			return false;
		}
		inline bool testRecvExtraTag(int src, int* tag)
		{
			int flag;
			MPI_Status status;
			MPI_Test(&recvRequests[src], &flag, &status);
			if(flag){
				*tag = status.MPI_TAG;
				return true;
			}
			return false;
		}

		inline bool testSend(int dest){
			int flag;
			MPI_Test(&sendRequests[dest], &flag, MPI_STATUS_IGNORE);
			if(flag)
				return true;
			else
				return false;
		}
		inline void cancelSend(const int dest){
			MPI_Cancel(&sendRequests[dest]);
		}
		inline void cancelRecv(const int src){
			MPI_Cancel(&recvRequests[src]);
		}
		inline void cancelRecv(){
			MPI_Cancel(&recvRequests[_rank]);
		}

		inline bool iprobe(const int src = MPI_ANY_SOURCE, const int tag = MPI_ANY_TAG)
		{
			int flag;
			MPI_Iprobe(src, tag, comm, &flag, MPI_STATUS_IGNORE);
			if(flag)
				return true;
			else
				return false;
		}
		inline void barrier()
		{
			MPI_Barrier(comm);
		}
		//
		//Advanced Operations
		//
		void allReduce(void* sendBuf, void* recvBuf, int count, DataType_t datatype, const Op_t op){
			MPI_Allreduce(sendBuf, recvBuf, 1, datatype, op, comm);
		}
		int reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root){
			return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
		}
		//
		inline int voteToHalt(const bool isActive=false){
			int ret, myState;
			myState = isActive? 1:0;
			allReduce(&myState, &ret, 1, MPI_INT, MPI_SUM);
			return ret;
		}
		//
		template <typename T>
		void allGather(T& elem, std::vector<T>& results){
			//Not implemented.
			//MPI_Allgather(sendBuf, sendCount, datatype, recvBuf, recvCount, comm);
		}
		//
		template <typename T>
		void gather(T& elem, std::vector<T>& results){
			//Not implemented.
		}
		//timer
		inline void startTimer(){
			time = MPI_Wtime();
		}
		inline double stopTimer(){
			const double tmp = time;
			time = MPI_Wtime();
			return (time - tmp);
		}
		//profiling
		void initProfiling(){
			#ifdef Comm_Profiling_Send
			sendDataVolume = 0;
			#endif
			#ifdef Comm_Profiling_Recv
			recvDataVolume = 0;
			#endif
		}
		void outputProfilingData(const int root = 0){
			//barrier();
			#ifdef Comm_Profiling_Send
			long vol1 = 0;
			reduce(&sendDataVolume, &vol1, 1, MPI_LONG, MPI_SUM, root);
			if(_rank == root){
				std::cout<<"Totoal Send Data(Byte): "<<vol1<<std::endl;
			}
			#endif
			#ifdef Comm_Profiling_Recv
			long vol2 = 0;
			reduce(&recvDataVolume, &vol2, 1, MPI_LONG, MPI_SUM, 0);
			if(rank == root){
				std::cout<<"Totoal Recv Data(Byte): "<<vol2<<std::endl;
			}
			#endif
		}
	};
	}//end COMM
}
#endif
