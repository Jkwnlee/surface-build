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

def remove_overlap(list_a):
	temp=[]
	list_a = sorted(list_a)
	for a in range(len(list_a)):
		if a >= ( len(list_a) - 2 ):
			pass
		elif JH_lib.distance_atoms(list_a[a] ,list_a[a+1] ) > perterbation : 
			temp.append(list_a[a])
	return temp

def generate_add_ads(site, range_from_adsorp, orig_slab, outmost_slab, orig_unitcell, orig_adsorp, name_of_site):
	i=1
	temp_cell, temp_compound, expand_slab = JH_lib.pst_cell_expansion(unitcell, outmost_slab, [4,4,1])
	expand_slab = JH_lib.pst_poscar_move(expand_slab, float(unitcell[0][0])+float(unitcell[1][0]), float(unitcell[1][1]), 0)

	for a in site: 
		## initialization of list "adsorp_position"
		adsorp_position=[]
		for b in orig_adsorp: adsorp_position.append(b)
		## ----------------------------------------------

		if range_from_adsorp[0] < JH_lib.distance_atoms(orig_adsorp[0][1:4], a) <  range_from_adsorp[1] :
			## Calculate normal vector and shift with following the vector
			near_atoms=[]; j = 1; temp_vec=[]; moved_site=[]
			if name_of_site == 'top':
				for b in range(len(a)):
					# print b
					if b != 3:	moved_site.append(a[b])
					else:	moved_site.append(a[b]+range_from_adsorp[0]*1.5)
			else:
				for b in range(len(a)): moved_site.append(a[b])
				for normal in expand_slab:
					if JH_lib.distance_atoms(a, normal[1:4]) <  range_from_adsorp[1] and j < 4:
						# near_atoms.append([normal[1],normal[2],expand_slab[0][2]])
						near_atoms.append(normal[1:4])
						j = j + 1
				normal_v, center_of_plane = JH_lib.normal_vector(near_atoms)
				# # for xxx in  near_atoms: print xxx
				# # print normal_v
				temp_vec, moved_site = JH_lib.shift_ads_w_normal_vec(a,range_from_adsorp[0], normal_v)
			## ----------------------------------------------

			print name_of_site, orig_adsorp[0][1:4], moved_site, JH_lib.distance_atoms(orig_adsorp[0][1:4], moved_site),range_from_adsorp
			adsorp_position.append(['C', moved_site[0],moved_site[1], moved_site[2], 'T', 'T', 'T'])
			JH_lib.add_adsorbate(orig_slab, orig_unitcell, adsorp_position, str(name_of_site)+'_'+str(i)+'.vasp')
			i = i + 1


## Test for normal vector 
# temp, center=normal_vector([[6.0598862556776316,0.0000000000000000,4.3528408899082995],[ 6.0042276279323659,2.5098483562000000,4.2973650164883086],[4.2408565966276779,1.2419402764485341,5.5002446586139300]])
# shift_ads_w_normal_vec(center,1.09601,temp)



# perterbation : distance mapping criterial (collecting the distances less than perterbation,  )
unitcell, compound, position = JH_lib.r_cryst_vasp('CONTCAR')
 
perterbation=1.5;

slab=[];adsorp=[];
new_distance=[]; distance=[]
outmost_al=[] ; N_bb=[]
sub_al=[]

temp=0; i = 1; j = 1; temp_high=0
for a in position:
	if a[0] == 'Cu':
		a[0]=a[0]+'_'+str(i)
		slab.append(a)
		i = i + 1
		# find highest surface atom
		if a[3] >  temp_high: temp_high = a[3]
	else:
		adsorp.append(a)

a1=[]; a2=[]
for a in slab:
	for b in slab:
		if a[1:4] != b[1:4]: distance.append([JH_lib.distance_atoms(a[1:4], b[1:4]), a[0], b[0]])
distance= sorted(distance)


'''
Collect the outmost layer with considering 
(1) the number of broken bonds (less than ideal)
and 
(2) higher than middle of slab
'''
dimension = [3, 3 , 1]
temp_cell, temp_compound, expand_slab = JH_lib.pst_cell_expansion(unitcell, slab, dimension)
expand_slab = JH_lib.pst_poscar_move(expand_slab, float(unitcell[0][0])+float(unitcell[1][0]), float(unitcell[1][1]), 0)
for a in slab:
	temp = 0
	for b in expand_slab:
		if JH_lib.distance_atoms(a[1:4], b[1:4]) <  distance[0][0] * 1.1 and a[0] != b[0]:
			temp = temp + 1
	N_bb.append(temp)


for a in range(len(N_bb)):
	if N_bb[a] < 12 and slab[a][3] > temp_high/2 : 
		outmost_al.append(slab[a])
	# if a[3] > ( temp_high - perterbation - distance[0][0]): sub_al.append(a)

## ------consider the sub surface later-----


'''
Survay possible site
'''


top=[]; temp_t=[]
bridge=[]; temp_b=[]
hollow=[]; temp_h=[]
hollow_4=[]; temp_h4=[]



for a in outmost_al:
	print a
	temp=[];	temp=a[1:4]
	temp[2] = temp[2] +perterbation
	temp_t.append(temp)

	# find top site (1 hand)
	for b in outmost_al:
		# find bridge site (2 hand)
		if a != b and JH_lib.distance_atoms(a[1:4],b[1:4]) <  (distance[0][0] + perterbation):
			temp_p1= a[1:4]
			temp_p2= b[1:4]
			temp=[]; temp=average([temp_p1,temp_p2])
			temp_b.append(temp)
			# find hollow site (3 bond)
			for c in outmost_al:
				if c != a and c != b and \
					JH_lib.distance_atoms(a[1:4],c[1:4]) <  (distance[0][0] + perterbation) and\
					JH_lib.distance_atoms(b[1:4],c[1:4]) <  (distance[0][0] + perterbation):
					temp_p3=c[1:4]
					temp=[];temp=average([temp_p1,temp_p2,temp_p3])
					temp_h.append(temp)
					# find 4 fold hollow site (4 bond)
					for d in outmost_al:
						if d != a and d != b and d != c and \
							JH_lib.distance_atoms(a[1:4],d[1:4]) <  (distance[0][0] + perterbation) and\
							JH_lib.distance_atoms(b[1:4],d[1:4]) <  (distance[0][0] + perterbation) and\
							JH_lib.distance_atoms(c[1:4],d[1:4]) <  (distance[0][0] + perterbation):
							temp_p4=d[1:4]
							temp=[];temp=average([temp_p1,temp_p2,temp_p3,temp_p4])
							temp_h4.append(temp)

top=remove_overlap(temp_t)
bridge=remove_overlap(temp_b)
hollow=remove_overlap(temp_h)
hollow_4=remove_overlap(temp_h4)


print hollow

for a in range(len(slab)):
	slab[a][0]=  slab[a][0].split('_')[0]

generate_add_ads(top,      [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'top')
generate_add_ads(bridge,   [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'bridge')
generate_add_ads(hollow,   [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'hollow')
generate_add_ads(hollow_4, [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'hollow4f')
