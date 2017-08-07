#!/bin/sh
coverage=$1
system=$2
b=`echo $coverage"_"$system`

rm output.txt
cd $b
for a in *
do cd $a
echo $b $a `tail -1 OSZICAR | awk '{printf "%10.9f", $5}'` >> ../../output.txt
cd ..
done
cd ..

