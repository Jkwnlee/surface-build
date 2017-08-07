#!/bin/python
import os
import math 
import sys
import time
import JH_lib
import pylab

def average(input_vectors):
	x=0; y=0; z=0
	if len(input_vectors[0]) == 3:
		for a in input_vectors:
			x = x + float(a[0]) ; y = y + float(a[1]); z = z + float(a[2])

		return [x/len(input_vectors), y/len(input_vectors) ,z/len(input_vectors) ]
	else:
		JH_lib.print_error(True,'average')

def remove_overlap(list_x, unitcell):
	list_a=[]
	for a in list_x:
		if a[0] <  ( float(unitcell[0][0]) + float(unitcell[1][0]) + float(unitcell[2][0])) and\
		 a[1] <  ( float(unitcell[0][1]) + float(unitcell[1][1]) + float(unitcell[2][1])) and\
		 a[2] <  ( float(unitcell[0][2]) + float(unitcell[1][2]) + float(unitcell[2][2])):
			list_a.append(a)

	list_a = sorted(list_a); 	temp=[list_a[0]]; 	i = 0 
	for a in range(len(list_a)):
		for b in range(len(list_a)):
			if a == b or list_a[a] == list_a[b] or list_a[a] == temp[i]: pass
			else: temp.append(list_a[a]); i = i + 1
	return sorted(temp)


def generate_add_ads(site, range_from_adsorp, center, orig_slab, orig_unitcell, orig_adsorp, name_of_site):
	i=1

	for a in site: 
		center[2] = a[2]
		if range_from_adsorp[0] < JH_lib.distance_atoms(a, center) < range_from_adsorp[1]:
		## initialization of list "adsorp_position"
			adsorp_position=[]; moved_site=[]
			moved_site = JH_lib.pst_poscar_move(orig_adsorp,  - a[0], - a[1], 0)
			JH_lib.add_adsorbate(orig_slab, orig_unitcell, moved_site, str(name_of_site)+'_'+str(i)+'.vasp')
			moved_site = JH_lib.pst_poscar_move(orig_adsorp,  a[0], a[1], 0)

			i = i + 1

