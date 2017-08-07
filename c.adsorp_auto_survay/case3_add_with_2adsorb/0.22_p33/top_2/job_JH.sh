#!/bin/sh
#PBS -q sandy
#PBS -l nodes=4:ppn=12
#PBS -l walltime=999:00:00

NPROCS=`wc -l < $PBS_NODEFILE`

echo "PBS_O_WORKDIR = ${PBS_O_WORKDIR}"
echo "PBS_JOBID     = ${PBS_JOBID}"

hostname
date

cd $PBS_O_WORKDIR
cp $PBS_NODEFILE nodefile

function run_vasp(){
module add Intel/Compiler/17.0.1
module add Intel/MKL/2017.1-132
module add Intel/MPI/2017.1-132
## Please check the binary files in /GRAPE/Apps/VASP/OLD_BIN/
## vasp 5.3.5 version
mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/OLD_2015_RHEL6_BIN/vasp.5.3.5_31MAR2014_GRAPE_GRAPHENE_NORMAL.mpi.x > stdout.log
}


run_vasp

exit 0
