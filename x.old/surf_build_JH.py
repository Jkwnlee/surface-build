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

###########################################################################
###########################################################################
####   ____   ___    _   _  ___ _____   _____ ___  _   _  ____ _   _   ####
####  |  _ \ / _ \  | \ | |/ _ \_   _| |_   _/ _ \| | | |/ ___| | | |  ####
####  | | | | | | | |  \| | | | || |     | || | | | | | | |   | |_| |  ####
####  | |_| | |_| | | |\  | |_| || |     | || |_| | |_| | |___|  _  |  ####
####  |____/ \___/  |_| \_|\___/ |_|     |_| \___/ \___/ \____|_| |_|  ####
####                                                                   ####
###########################################################################
###########################################################################
def readPOSCAR(filename):
    """
    Filename: VASP, POSCAR TYPE
    Function: Read POSCAR type 
    """
    poscar=[]; unitcell=[]; compound=[]; system=[];scales=[]
    with open(filename, 'r') as f:
        i = 1; j = 1; k =0
        for line in f:
            if i > 2 and i <6:
                unitcell.append(line.split())
            elif i > 5 and i < 8:
                compound.append(line.split())
            elif i == 8 :
                if j == 1 :
                    if line.split()[0][0] == 'S':
                        # print 'Selective Dynamics are applied'
                        i = i - 1
                    j = 2
                if j == 2:
                    if line.split()[0][0] == 'C':
                        scales=[1.,1.,1.]
                    elif line.split()[0][0] == 'D':
                        scales=[float(unitcell[0][0]) , float(unitcell[1][1]), float(unitcell[2][2])]
            elif i > 8 :
                x = scales[0]  * float(line.split()[0])
                y = scales[1]  * float(line.split()[1])
                z = scales[2]  * float(line.split()[2])
                if k < int(compound[1][0]):
                    system.append([compound[0][0], x, y, z])
                elif int(compound[1][0]) <= k <= int(compound[1][1]):
                    system.append([compound[0][1], x, y, z])
                k= k+1
            i = i + 1
    poscar=(unitcell, compound, system)
    return poscar, unitcell, compound, system

def poscar_head(num_atom,num_rlx):
    scale=( num_atom - 1 ) / 2
    num_fix_lay=( num_atom - num_rlx ) / 2 


def judge_pri_conv(filename, system):
    ## Function check the condition of file wheather it is primitive or not
    ##
    ## Condition1: number atom == 1
    ## Condition2: nit
    cell, unitcell, compound, systemll=readPOSCAR('CONTCAR')
    if len(compound[1]) == 1:
        # print "this cell include only one kind of atom"
        if compound[1][0] == '1':
            # print "this cell include only one of atom"
            True
    else:
        print "this cell include only more than one kind of atoms", compound[1]
        pass

    print unitcell
    if system == 'fcc':
        material_name=compound[0][0]
        latcon=unitcell[0][0]
        if float(latcon) == 0.:
            # print "the unit cell is not conventional"
            latcon=unitcell[0][1] * 2
        else:
            pass
    return latcon


def write_poscar(poscar, filename):
    """
    write a poscar file type
    """

# mtg=mtg_structure

cell, unitcell, compound, systemll=readPOSCAR('CONTCAR')
system = 'fcc'
judge_pri_conv(filename, system):

write_poscar('test', cell)

print compound[1]




###########################################################################
###########################################################################
####   ____   ___    _   _  ___ _____   _____ ___  _   _  ____ _   _   ####
####  |  _ \ / _ \  | \ | |/ _ \_   _| |_   _/ _ \| | | |/ ___| | | |  ####
####  | | | | | | | |  \| | | | || |     | || | | | | | | |   | |_| |  ####
####  | |_| | |_| | | |\  | |_| || |     | || |_| | |_| | |___|  _  |  ####
####  |____/ \___/  |_| \_|\___/ |_|     |_| \___/ \___/ \____|_| |_|  ####
####                                                                   ####
###########################################################################
###########################################################################
# target_crystal='fcc' #'hcp'

# target_material='Pt'    #'Os'
# a_len=2.7323459964151184
# # c_len=4.3241327206235365
# c_len=a_len
# hkl_index=(1, 1, 3)
# vacuum_distance=18
# number_of_layer = 181  ## Unit: Atomic layer or number of atom

# ## select_dynamic[0] = True (Apply)
# ##                   = False (Skip this option)
# ## select_dynamic[1] = Thickness of bulk-like region
# ## select_dynamic[2] = Thickness of slab
# select_dynamic=['True', 6 , 15 ]  
# re_arrange='True'  # = True (Apply) ; False (Skip this option)
# filename=target_material+'_surf'


