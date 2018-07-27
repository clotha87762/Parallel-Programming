
# NOTICE: Please do not remove the '#' before 'PBS'

# Name of your job
#PBS -N MPI_dynamic_weak

# Declaring job as not re-runnable
#PBS -r n

# Resource allocation (how many nodes? how many processes per node?)
#PBS -l nodes=4:ppn=12

#PBS -o outputDynamic2.txt
#PBS -e errorDynamic2.txt
# Max execution time of your job (hh:mm:ss)
# Your job may got killed if you exceed this limit
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
NUM_MPI_PROCESS_PER_NODE=4	# edit this line to set number of MPI process you want to use per node
export OMP_NUM_THREADS=3	# set max number of threads OpenMP can use per MPI task
export MV2_ENABLE_AFFINITY=0	# prevent MPI from binding all threads to one core
# Please make sure that NUM_MPI_PROCESS_PER_NODE * OMP_NUM_THREADS == ppn

time mpiexec -ppn $NUM_MPI_PROCESS_PER_NODE ./MS_Hybrid_dynamic2 $OMP_NUM_THREADS -2 2 -2 2 2000 2000 disable # edit this line to fit your needs!
# In this case, it will use 3 nodes to run, each nodes with 3 processes, each process with 4 threads. (3 process * 4 thread = 12 ppn )
# Besides NUM_MPI_PROCESS_PER_NODE, OMP_NUM_THREADS, you also need to modify the number of nodes or ppn you want to use.
# 7 10 14 20 28 40 57 80 113 160 226 