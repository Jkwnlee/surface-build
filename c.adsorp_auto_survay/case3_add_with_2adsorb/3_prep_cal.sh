#!/bin/sh

## Ji-Hwan's functions
shell_path="/home/jihwan/code_work/shell_script"
. $shell_path/kohn_rewrite_job.sh


system=$1
low_system=$2

cp $system/$low_system/CONTCAR .

cd vasp

for structure in *.vasp;do

temp=$(echo $structure | cut -c 1-`echo ${#structure} - 5 | bc` )
mkdir $temp; mv $structure $temp/POSCAR
cp ../$system/$low_system/{POTCAR,INCAR,KPOINTS,job*sh} $temp/.

cd $temp
rewrite_job -g job*.sh sandy 4; mv new* job*.sh;
#qsub job*.sh
cd ..

done

cd ..
