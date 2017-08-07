#!/bin/python

### The python code for atomatic Wulff construction
### Written by Ji-Hwan Lee
### Last update:20 jan 2017
### Used library:os, math, ase
### Purpose: make wulff construction  --> go to 
### step 1: Revise surface_slab_fcc_bcc into python
### setp 2: 

### The python code for 
### Written by Ji-Hwan Lee
### Please use this code and cite our surface energy library
### S.-H. Yoo, J.-H. Lee, Y.-K. Jung, and A. Soon, Phys. Rev. B 93, 035434 (2016).

import os
import math 
import sys
import time


## ------------------------------------------------------------------------------
"""
LIST-A [FOR STRUCTURES]
1. r_cryst_vasp
: reading VASP crystal structure file type 
(.vasp or called 'POSCAR, CONTCAR')
"""
## ------------------------------------------------------------------------------


def r_cryst_vasp(filename):
    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    poscar=[]; unitcell=[]; compound=[]; position=[];scales=[];num_atoms=0
    with open(filename, 'r') as f:
        i = 1; j = 1; k =0; Selective=False; kn = 0
        for line in f:
            if line == None: break
            if i > 2 and i <6:
                unitcell.append(line.split())
            elif i > 5 and i < 8:
                compound.append(line.split())
            elif i == 8 :
                if j == 1 :
                    if line.split()[0][0] == 'S':
                        # print 'Selective Dynamics are applied'
                        Selective= True
                        i = i - 1
                    j = 2
                if j == 2:
                    if line.split()[0][0] == 'C':
                        scales=[[1.,0,0],[0,1.,0],[0,0,1.]]
                    elif line.split()[0][0] == 'D':
                        scales=[[ float(unitcell[0][0]), float(unitcell[0][1]), float(unitcell[0][2])],\
                                [ float(unitcell[1][0]), float(unitcell[1][1]), float(unitcell[1][2])],\
                                [ float(unitcell[2][0]), float(unitcell[2][1]), float(unitcell[2][2])]]
                if num_atoms == 0:
                    for temp in compound[1]: num_atoms=num_atoms+int(temp)

            elif i > 8 :
                if i <= 9 + num_atoms:
                    x = scales[0][0] * float(line.split()[0]) + scales[1][0] * float(line.split()[1]) + scales[2][0] * float(line.split()[2])
                    y = scales[0][1] * float(line.split()[0]) + scales[1][1] * float(line.split()[1]) + scales[2][1] * float(line.split()[2])
                    z = scales[0][2] * float(line.split()[0]) + scales[1][2] * float(line.split()[1]) + scales[2][2] * float(line.split()[2])
                    if k <= int(compound[1][kn]):
                        if Selective:
                            position.append([compound[0][kn], x, y, z, line.split()[3], line.split()[4], line.split()[5]])
                        else:
                            position.append([compound[0][kn], x, y, z])
                    k= k+1

                    if k == int(compound[1][kn]):
                        kn = kn + 1
                        k = 0


            if i == 8 + num_atoms:
                return unitcell, compound, position
            else:
                i = i + 1


## ------------------------------------------------------------------------------


def counts(a, b, c, x, y, z, length):
    if x * length <= a < (x + 1) * length and \
       y * length <= b < (y + 1) * length and \
       z * length <= c < (z + 1) * length: 
        return 1.
    else:
    #     # print a, b, c, x, y, z,0
        return 0.

## TEST for pst_cell_expansion
unitcell, compound, position=r_cryst_vasp('t300K_POSCAR')

Volume = float(unitcell[0][0]) * float(unitcell[1][1]) * float(unitcell[2][2])

scale=3
xwmax = int(float(unitcell[0][0]) /scale + 1)
ywmax = int(float(unitcell[1][1]) /scale + 1)
zwmax = int(float(unitcell[2][2]) /scale + 1)
grid = []; tarbulate = []
for zw in range(zwmax):
    for yw in range(ywmax):
        for xw in range(xwmax):
            grid.append([float(xw), float(yw), float(zw)])
            tarbulate.append(float(0))

system = []
for a in position:
    if a[0] == 'O':
        if a[1] < 0:
            x = a[1] + float(unitcell[0][0])
        else:
            x = a[1]
        if a[2] < 0:
            y = a[2] + float(unitcell[1][1])
        else:
            y = a[2]
        if a[2] < 0:
            z = a[3] + float(unitcell[2][2])
        else:
            z = a[3]
        system.append(['O', x, y, z])

for a in system:
    for b in range(len(grid)):
        tarbulate[b] = tarbulate[b] + counts(a[1], a[2],a[3], grid[b][0], grid[b][1], grid[b][2], scale)

i = 0
print xwmax, ywmax, zwmax
for a in tarbulate:
    print "%18.10e" %( a * Volume ) ,
    if i == 4:
        i = 0
        print ""
    else:
        i = i + 1









