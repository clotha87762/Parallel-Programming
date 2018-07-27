#PBS -N My_JOB
#PBS -r n
#PBS -l nodes=1:ppn=5
#PBS -l walltime=00:05:00
#PBS -e errorADV.txt
#PBS -o outputADV.txt

cd $PBS_O_WORKDIR
mpiexec ./advanced 15 myinput myAdvOutput
