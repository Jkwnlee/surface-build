#!/bin/sh
system=$1
low_system=$2

cp $system/$low_system/CONTCAR .

echo "Choose the most stable system, "$low_system
python ./auto_survay.py
