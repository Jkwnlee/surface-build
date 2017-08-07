#!/bin/python
import os
import math 
import sys
import time
import subprocess



coverage = 0.25
periodic = 'p22'
print  '%10.5s %10.3s' %(str(coverage), periodic)
os.system('./1_output_e.sh %10.5s %10.5s' %(str(coverage), periodic) )

output=[]; low_system=[]
with open('output.txt', 'r') as f:
	for a in f: 
		output.append(a.split())

for a in output:
	for b in output:
		if a == b:
			pass
		elif float(a[2]) < float(b[2]):
			low_system=a

os.system('./2_copy_contcar.sh %10.9s %10.5s' \
	%(low_system[0], low_system[1]) )

