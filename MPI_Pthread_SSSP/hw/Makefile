CC			:= icc
CXX			:= icpc

MPICC			:= mpiicc
MPICXX			:= mpiicpc

CCFLAGS			:= -O3 -xHost -Wall -std=gnu11
CXXFLAGS		:= -O3 -xHost -Wall -std=c++0x
LDFLAGS			:= -lX11 -fopenmp -pthread 

all: SSSP_Pthread SSSP_MPI_sync SSSP_MPI_async 


SSSP_Pthread: SSSP_Pthread.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?	
SSSP_MPI_sync: SSSP_MPI_sync.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
SSSP_MPI_async: SSSP_MPI_async.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
SSSP_Pthread_test: SSSP_Pthread_test.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?	
SSSP_MPI_sync_test: SSSP_MPI_sync_test.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
SSSP_MPI_async_test: SSSP_MPI_async_test.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
clean:
	rm -f SSSP_Pthread SSSP_MPI_sync SSSP_MPI_async SSSP_Pthread_test SSSP_MPI_sync_test SSSP_MPI_async_test
