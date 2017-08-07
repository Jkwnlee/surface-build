#!/bin/python
import os
import math 
import sys
import time
import JH_lib

## Example for rotation_mlc to generate  [OUTPUT1]
# for angle in range(10):
#     temp = 180./float(angle+1)
#     JH_lib.rotate_mlc('grycerol.vasp', temp, 'z', True)

## Example for adsorp it on the surface [OUTPUT2]
# slab_unitcell, compound, slab_position = JH_lib.r_cryst_vasp('slab_rh44.vasp') # READ Slab
# adsorbate_position = JH_lib.mlc_grycerol([0.0000000000000000,0.0000000000000000,10.2074228061374646])
# JH_lib.add_adsorbate(slab_position, slab_unitcell, adsorbate_position, 'output2.vasp')


# Example for adsorp/rotate on surface [OUTPUT3]
# position_of_ref=[0.0000000000000000,0.0000000000000000,11.2074228061374646]
# slab_unitcell, compound, slab_position = JH_lib.r_cryst_vasp('slab_rh44.vasp') # READ Slab

# for angle in range(10):
# 	adsorbate_position = JH_lib.mlc_grycerol(position_of_ref)
# 	temp = 180./float(angle+1)
# 	compound, new_position = JH_lib.rotate_mlc(adsorbate_position,  temp, 'z')
# 	JH_lib.add_adsorbate(slab_position, slab_unitcell, new_position, "output3_"+str(angle)+'.vasp')


