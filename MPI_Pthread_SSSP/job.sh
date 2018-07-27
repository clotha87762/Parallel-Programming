
# NOTICE: Please do not remove the '#' before 'PBS'

# Name of your job
#PBS -N HYBRID

# Declaring job as not re-runnable
#PBS -r n

# Resource allocation (how many nodes? how many processes per node?)
#PBS -l nodes=1:ppn=12

#PBS -o print.txt
#PBS -e MPI_async500/error1_12.txt

#PBS -l walltime=00:03:00

cd $PBS_O_WORKDIR

export MV2_ENABLE_AFFINITY=0	# prevent MPI from binding all threads to one core

time mpirun -np 500 ./SSSP_MPI_async_test 1 input500.txt output.txt 250 # edit this line to fit your needs!
