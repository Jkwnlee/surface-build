#!/bin/sh
#PBS -q west2
#PBS -l nodes=4:ppn=12
####PBS -l walltime=48:00:00

NPROCS=`wc -l < $PBS_NODEFILE`

echo "PBS_O_WORKDIR = ${PBS_O_WORKDIR}"
echo "PBS_JOBID     = ${PBS_JOBID}"

hostname
date

cd $PBS_O_WORKDIR
cp $PBS_NODEFILE nodefile


bin_path=/home/emin/programs/lammps-17Nov16/src/lmp_mpi

function run_lammps(){
module add Intel/Compiler/17.0.1
module add Intel/MKL/2017.1-132
module add Intel/MPI/2017.1-132

mpirun -genv I_MPI_DEBUG 5 -np $NPROCS $bin_path < in.min

}
for a in c*
do cd $a
run_lammps
cd ..
done
exit 0