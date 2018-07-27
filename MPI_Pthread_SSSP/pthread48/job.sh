
# NOTICE: Please do not remove the '#' before 'PBS'

# Name of your job
#PBS -N HYBRID

# Declaring job as not re-runnable
#PBS -r n

# Resource allocation (how many nodes? how many processes per node?)
#PBS -l nodes=4:ppn=12

#PBS -o print.txt
#PBS -e pthread_30000/error100.txt
# Max execution time of your job (hh:mm:ss)
# Your job may got killed if you exceed this limit
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
NUM_MPI_PROCESS_PER_NODE=3	# edit this line to set number of MPI process you want to use per node
export NUM_THREADS=100	# set max number of threads OpenMP can use per MPI task
export MV2_ENABLE_AFFINITY=0	# prevent MPI from binding all threads to one core
# Please make sure that NUM_MPI_PROCESS_PER_NODE * OMP_NUM_THREADS == ppn

time mpiexec -np 48 ./SSSP_Pthread_test $NUM_THREADS intput30000sparse.txt output.txt 15000 # edit this line to fit your needs!
# In this case, it will use 3 nodes to run, each nodes with 3 processes, each process with 4 threads. (3 process * 4 thread = 12 ppn )
# Besides NUM_MPI_PROCESS_PER_NODE, OMP_NUM_THREADS, you also need to modify the number of nodes or ppn you want to use.
