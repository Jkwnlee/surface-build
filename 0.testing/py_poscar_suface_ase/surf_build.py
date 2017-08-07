#!/bin/python

### The python code for atomatic Wulff construction
### Written by Ji-Hwan Lee
### Last update:20 jan 2017
### Used library:os, math, ase
### Purpose: make wulff construction  --> go to 
### step 1: 
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

def MTG_surf_buiding(target_crystal,
	target_material, a_len, c_len, \
	hkl_index, vacuum_distance, number_of_layer, \
	select_dynamic, re_arrange, filename):

	from MTG_post_ase import post_ase_vasp_rewrite

	from ase.io import vasp, write
	from ase.lattice.cubic import FaceCenteredCubic
	from ase.lattice.hexagonal import HexagonalClosedPacked
	from ase import Atoms, Atom
	from ase.lattice.surface import surface
	# from general_surface import surface
	# from ase.lattice import general_surface


	if target_crystal == 'hcp':
		hcp_bulk = HexagonalClosedPacked(symbol=target_material, 
			latticeconstant={'a':a_len, 'c':c_len}
			)  
		hcp_surf = surface(hcp_bulk, hkl_index , number_of_layer/2)
		hcp_surf.center(vacuum=vacuum_distance/2, axis=2)
		vasp.write_vasp(filename+'_temp.vasp', hcp_surf)
		post_ase_vasp_rewrite(filename+'_temp.vasp',filename+'.vasp', re_arrange, vacuum_distance, select_dynamic)#,'False')

	elif target_crystal == 'fcc':
		fcc_bulk = FaceCenteredCubic(
			symbol            = target_material, 
			latticeconstant   = a_len
			)  
		if number_of_layer % 2. == 0 :
			number_of_layer= number_of_layer+1
			print """the input number_of_layer is even number.\n \
For FCC system, we recommend to use odd number of number_of_layer. \n \
Automatically 1 more layer is added, so the number_of_layer is """,number_of_layer

		else:
			print number_of_layer /2.
		fcc_surf = surface(fcc_bulk, hkl_index , number_of_layer)#, vacuum=vacuum_distance)
		fcc_surf.center(vacuum=vacuum_distance/2, axis=2)
		vasp.write_vasp(filename+'_temp.vasp', fcc_surf)
		post_ase_vasp_rewrite(filename+'_temp.vasp',filename+'.vasp', re_arrange, vacuum_distance, select_dynamic)#,'False')

	


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

target_crystal='fcc' #'hcp'

target_material='Pt'    #'Os'
a_len=2.7323459964151184
# c_len=4.3241327206235365
c_len=a_len
hkl_index=(1, 1, 3)
vacuum_distance=18
number_of_layer = 181  ## Unit: Atomic layer or number of atom

## select_dynamic[0] = True (Apply)
##                   = False (Skip this option)
## select_dynamic[1] = Thickness of bulk-like region
## select_dynamic[2] = Thickness of slab
select_dynamic=['True', 6 , 15 ]  
re_arrange='True'  # = True (Apply) ; False (Skip this option)
filename=target_material+'_surf'



MTG_surf_buiding(target_crystal,
	target_material, a_len, c_len, \
	hkl_index, vacuum_distance, number_of_layer, \
	select_dynamic, re_arrange, filename)
