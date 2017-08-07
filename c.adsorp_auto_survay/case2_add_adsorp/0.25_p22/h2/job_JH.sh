#!/bin/sh
#PBS -q west
#PBS -l nodes=4:ppn=8
#PBS -l walltime=48:00:00

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

#---------------- vasp 5.3.5 ----------------#
mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.3.5/NORMAL/vasp.5.3.5_31MAR2014_GRP7_NORMAL_VTST.x > stdout # For other nodes
#mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.3.5/AVX1/vasp.5.3.5_31MAR2014_GRP7_AVX1.x > stdout # For sandy nodes
#mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.3.5/AVX2/vasp.5.3.5_31MAR2014_GRP7_AVX2_VTST.x > stdout # For xeon nodes

#---------------- vasp 5.4.1 ----------------#
#mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.4.1/NORMAL/vasp_5.4.1_GRP7_NORMAL_p13082016.x # For other nodes
#mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.4.1/AVX1/vasp_5.4.1_GRP7_AVX1_p13082016.x > stdout # For sandy nodes
#mpirun -genv I_MPI_DEBUG 5 -np $NPROCS /GRAPE/Apps/VASP/bin/5.4.1/AVX2/vasp_5.4.1_GRP7_AVX2_p13082016.x > stdout # For xeon nodes


}


run_vasp

exit 0
