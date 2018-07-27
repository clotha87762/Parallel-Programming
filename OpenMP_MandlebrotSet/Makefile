CC			:= gcc
CXX			:= icpc

MPICC			:= mpicc
MPICXX			:= mpiicpc

CCFLAGS			:= -O3 -march=native -Wall -std=gnu11
CXXFLAGS		:= -O3 -xHost -Wall -std=c++0x
LDFLAGS			:= -lX11 -fopenmp -lm 

all: MS_MPI_static MS_MPI_dynamic MS_OpenMP_static MS_OpenMP_dynamic MS_Hybrid_static MS_Hybrid_dynamic MS_Hybrid_dynamic2 MS_seq

MS_MPI_static: MS_MPI_static.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_MPI_dynamic: MS_MPI_dynamic.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_OpenMP_static: MS_OpenMP_static.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_OpenMP_dynamic: MS_OpenMP_dynamic.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_Hybrid_static: MS_Hybrid_static.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_Hybrid_dynamic: MS_Hybrid_dynamic.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_Hybrid_dynamic2: MS_Hybrid_dynamic2.c
	$(MPICC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
MS_seq: MS_seq.c
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $?
clean:
	rm -f MS_MPI_static MS_MPI_dynamic MS_OpenMP_static MS_OpenMP_dynamic MS_Hybrid_static MS_Hybrid_dynamic MS_Hybrid_dynamic2 MS_seq 
