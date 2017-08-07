#!/bin/python
import os
import math 
import sys
import time
import subprocess



coverage = 0.22
periodic = 'p33'
print  'You choose the system r%10.5s %10.3s' %(str(coverage), periodic)
os.system('./1_output_e.sh %10.5s %10.5s' %(str(coverage), periodic) )

output=[]; low_system=[]
with open('output.txt', 'r') as f:
	for a in f: 
		output.append(a.split())

low_system=['system','cases',.999]
for a in range(len(output)):
	if float(output[a][2]) < float(low_system[2]):
		low_system =  output[a]

os.system('./2_copy_contcar.sh %10.9s %10.5s' \
	%(low_system[0], low_system[1]) )

