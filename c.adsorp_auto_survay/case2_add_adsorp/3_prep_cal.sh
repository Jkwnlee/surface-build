#!/bin/sh
system=$1
low_system=$2

cp $system/$low_system/CONTCAR .
python ./auto_survay.py
