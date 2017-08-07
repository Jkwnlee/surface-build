#!/bin/python
import os
import math 
import sys
import time
import JH_lib
import pylab
from  auto_survay import *
# -----------------------------------------------------------------------------
# Considering top/ bridge/ 4f-hollow site

unitcell, compound, position = JH_lib.r_cryst_vasp('input')
range_from_center = [0,3] 
# [0,3] From center, drowing circles with radious of 0 and 3 Angstrom where we search
distance = 3.2
v_site   = [6.012599744, 6.012599744,6.4486583697499995]
# -----------------------------------------------------------------------------
#                                 DO NOT TOUCH  B
#                                 DO NOT TOUCH  E
#                                 DO NOT TOUCH  L
#                                 DO NOT TOUCH  O
#                                 DO NOT TOUCH  W
# -----------------------------------------------------------------------------


slab=[];adsorp=[];
new_distance=[]; 
outmost_al=[] ; N_bb=[]
sub_al=[]

temp=0; i = 1; j = 1; temp_high=0; temp_low=float(unitcell[2][2])
for a in position:
	if a[0] == 'Ti' or a[0] == 'C':
		a[0]=a[0]+'_'+str(i)
		slab.append(a)
		i = i + 1
		# find highest surface atom
		if a[3] >  temp_high: temp_high = a[3]
	else:
		adsorp.append(a)
		# find lowest adatom
		if a[3] <  temp_low: temp_low = a[3]

perterbation = temp_low -  temp_high

adsorp = JH_lib.pst_poscar_move(adsorp,adsorp[1][1],adsorp[1][2],0)

'''
Collect the outmost layer with considering 
(1) the number of broken bonds (less than ideal)
and 
(2) higher than middle of slab
'''
dimension = [3, 3, 1]
temp_cell, temp_compound, expand_slab = JH_lib.pst_cell_expansion(unitcell, slab, dimension)
expand_slab = JH_lib.pst_poscar_move(expand_slab, float(unitcell[0][0])+float(unitcell[1][0]), float(unitcell[1][1]), 0)

for a in slab:
	if a [3] > 6:
		outmost_al.append(a)
## ------consider the sub surface later-----
'''
Survay possible site
'''

top=[]; temp_t=[]
bridge=[]; temp_b=[]
# hollow=[]; temp_h=[]
hollow_4=[]; temp_h4=[]

for a in outmost_al:
	# find top site (1 hand)
	temp=[];	temp=a[1:4]
	temp[2] = temp[2] + perterbation
	temp_t.append(temp)
	for b in outmost_al:
		# find bridge site (2 hand)
		if a != b and JH_lib.distance_atoms(a[1:4],b[1:4]) <  distance:
			temp_p1= a[1:4]
			temp_p2= b[1:4]
			temp=[]; temp=average([temp_p1,temp_p2])
			temp_b.append(temp)
			# find hollow site (3 bond)
			for c in outmost_al:
				if c != a and c != b and \
					JH_lib.distance_atoms(a[1:4],c[1:4]) <  distance and\
					JH_lib.distance_atoms(b[1:4],c[1:4]) <  distance:
					temp_p3=c[1:4]
					# temp=[];temp=average([temp_p1,temp_p2,temp_p3])
					# temp_h.append(temp)
					# find 4 fold hollow site (4 bond)
					for d in outmost_al:
						if d != a and d != b and d != c and \
							JH_lib.distance_atoms(a[1:4],d[1:4]) <   distance and\
							JH_lib.distance_atoms(b[1:4],d[1:4]) <   distance and\
							JH_lib.distance_atoms(c[1:4],d[1:4]) <   distance :
							temp_p4=d[1:4]
							temp=[];temp=average([temp_p1,temp_p2,temp_p3,temp_p4])
							temp_h4.append(temp)

top=remove_overlap(temp_t, unitcell)
bridge=remove_overlap(temp_b, unitcell)
# hollow=remove_overlap(temp_h, unitcell)
hollow_4=remove_overlap(temp_h4, unitcell)
# print len(top), len(temp_t)

generate_add_ads(top, range_from_center, v_site, slab, unitcell, adsorp, 'top')
generate_add_ads(bridge, range_from_center,  v_site, slab, unitcell, adsorp, 'bridge')
# generate_add_ads(hollow,   range_from_center,  v_site, slab, unitcell, adsorp, 'hollow')
generate_add_ads(hollow_4,  range_from_center,  v_site, slab, unitcell, adsorp, 'hollow4f')

os.system('mkdir vasp')
os.system('mv *.vasp vasp/.')
