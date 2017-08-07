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

def remove_overlap(list_a, unitcell, adsorp_position):
	i = 0 ; temp=[]
	list_a = sorted(list_a)
	temp=[list_a[0]]
	for a in range(len(list_a) - 1):
		if JH_lib.distance_atoms(list_a[a],list_a[a+1]) > 0.1:
			temp.append(list_a[a+1])
	print '		(*-1) Consider the ovelapped site , N_{site}: %3.i cases --> %3.i cases' %(len(list_a), len(temp))

	return temp


def bond_length_finding(input_position):
	a1=[]; a2=[]; distance=[]
	for a in range(len(input_position)):
		for b in range(a):
			if input_position[a][1:4] != input_position[b][1:4]: 
				distance.append(JH_lib.distance_atoms(input_position[a][1:4], input_position[b][1:4]))
	return min(distance)

def remove_overlap2_dist(adsorp, input_site):
	## Define the reference point 
	center_of_adsorp=[0,0,0]
	for a in adsorp:
		center_of_adsorp[0] = center_of_adsorp[0] + a[1]
		center_of_adsorp[1] = center_of_adsorp[1] + a[2]
		center_of_adsorp[2] = center_of_adsorp[2] + a[3]
	for a in range(3): center_of_adsorp[a] = center_of_adsorp[a]/len(adsorp)
	
	## removing sites which have same direction from the referecne point
	## [useful only for high - symmetric system]
	distance1=[];	distance2=[]
	criteria = [ float(unitcell[0][0]) + float(unitcell[1][0]) + float(unitcell[2][0]),\
				 float(unitcell[0][1]) + float(unitcell[1][1]) + float(unitcell[2][1]),\
				 0 ] 


 	for a in input_site:
		distance1.append([JH_lib.distance_atoms(a, center_of_adsorp),a])
 		temp2=[0,0,0]
		for b in range(3): temp2[b] = temp2[b] + center_of_adsorp[b] + criteria[b]
		distance2.append([JH_lib.distance_atoms(a, temp2),a])
	dist_from_ads1 = sorted(distance1)
	dist_from_ads2 = sorted(distance2)

	temp1=[dist_from_ads1[0][1]]
	for a in range(len(dist_from_ads1) - 1 ): 
		if ( dist_from_ads1[a + 1][0] - dist_from_ads1[a][0] ) > 0.1 : 
			temp1.append(dist_from_ads1[a+1][1])

	# dist_from_ads2 = sorted(distance2)
	# temp=[dist_from_ads2[0][1]]
	# for a in range(len(dist_from_ads) - 1 ): 
		# if ( dist_from_ads[a + 1][0] - dist_from_ads[a][0] ) > 0.1 : 
			# temp.append(dist_from_ads[a+1][1])
	print '''		(*-2) Skip atoms having same distance from adsorbate,
			  [useful only for high - symmetry system]
											N_{site}: %3.i cases --> %3.i cases''' %(len(input_site), len(temp1))#, len(temp2), len(temp)
	return temp1



def generate_add_ads(site, radious_from_adsroption, orig_position, outmost_slab, orig_unitcell, orig_adsorp, name_of_site):
	num_file = 1; range_from_adsorp=[0,radious_from_adsroption]
	temp_cell, temp_compound, expd_outmost_slab = JH_lib.pst_cell_expansion(orig_unitcell, outmost_slab, [3,3,1])
	expd_outmost_slab = JH_lib.pst_poscar_move(expd_outmost_slab, float(orig_unitcell[0][0])+float(orig_unitcell[1][0]), float(orig_unitcell[0][1]) +float(orig_unitcell[1][1]), 0)
	bond_criteria = bond_length_finding(outmost_slab) * 1.1
	dist = bond_length_finding(orig_position) * 0.9
	print '		[  ',
	for a in site: 
		## initialization of list "adsorp_position"
		adsorp_position=[]; i = 0
		while i < len(orig_adsorp):
			if dist < JH_lib.distance_atoms(a, orig_adsorp[i][1:4]) < radious_from_adsroption :	
				add_adsp=True
			else:	add_adsp=False
			i = i + 1

		if add_adsp:
			near_atoms=[]; j = 1; i = 0
			while j < 4:
				# print JH_lib.distance_atoms(expd_outmost_slab[i][1:4],a), bond_criteria, dist,bond_length_finding(outmost_slab)
				if JH_lib.distance_atoms(a, expd_outmost_slab[i][1:4]) <  bond_criteria:
					near_atoms.append(expd_outmost_slab[i][1:4])
					j = j + 1
				i = i + 1
			if len(near_atoms) != 3:
				return JH_lib.print_error(True,'the number of near_atoms is not 3 but %3.1i, please check code' %len(near_atoms))
			else:
				for b in range(3):
					for c in range(b):
						# print near_atoms[b][2],near_atoms[c][2],'before'
						if near_atoms[b][2] - near_atoms[c][2] > 1:
							pass
						else:
							near_atoms[c][2] = near_atoms[b][2]
						# print near_atoms[b][2],near_atoms[c][2],'after'

				normal_v, center_of_plane = JH_lib.normal_vector(near_atoms)
				temp_vec, moved_site = JH_lib.shift_ads_w_normal_vec(a,dist, normal_v)
				## ----------------------------------------------
				adsorp_position.append(['C', moved_site[0],moved_site[1], moved_site[2], 'T', 'T', 'T'])
				JH_lib.add_adsorbate(orig_position, orig_unitcell, adsorp_position, str(name_of_site)+'_'+str(num_file)+'.vasp') 
				print str(name_of_site)+'_'+str(num_file)+'.vasp ', 
				num_file = num_file + 1

	print '  ]  are generated'


# -----------------------------------------------------------------------------
filename='CONTCAR'
perterbation=1.5;
radious_from_adsroption=4.5
# -----------------------------------------------------------------------------
#                                 DO NOT TOUCH  B
#                                 DO NOT TOUCH  E
#                                 DO NOT TOUCH  L
#                                 DO NOT TOUCH  O
#                                 DO NOT TOUCH  W
# -----------------------------------------------------------------------------


slab=[];adsorp=[];
new_distance=[]; distance=[]
outmost_al=[] ; N_bb=[]
sub_al=[]

temp=0; i = 1; j = 1; temp_high=0
print '''# --------------------------- START --------------------------- #
This program is suppsed to useful only for the slab including %3.2s,
and adsorption of %3.2s
''' %('Cu','C')
unitcell, compound, position = JH_lib.r_cryst_vasp(filename)
print '(1) Reading input file, %10.9s' %filename

for a in position:
	if a[0] == 'Cu':
		a[0]=a[0]+'_'+str(i)
		slab.append(a)
		i = i + 1
		# find highest surface atom
		if a[3] >  temp_high: temp_high = a[3]
	else:
		adsorp.append(a)
bonds_length = bond_length_finding(slab)

print '''
(2) Reading the structure informations,
	(2-1) Seperate Slab part and Adsorption 
	(2-2) Check the highest atom in slab (%5.4f) in unit of Angstrom
	(2-3) Calculate the bond-length of slabs (%5.4f) in unit of Angstrom'''%(temp_high,bonds_length)
# -----------------------------------------------------------------------------



dimension = [3, 3, 1]
temp_cell, temp_compound, expd_slab = JH_lib.pst_cell_expansion(unitcell, slab, dimension)
expd_slab = JH_lib.pst_poscar_move(expd_slab, float(unitcell[0][0])+float(unitcell[1][0]), float(unitcell[1][1]), 0)
for a in slab:
	temp = 0
	for b in expd_slab:
		if JH_lib.distance_atoms(a[1:4], b[1:4]) <  bonds_length * 1.1 and a[0] != b[0]:
			temp = temp + 1
	N_bb.append(temp)

for a in range(len(N_bb)):
	if N_bb[a] < max(N_bb) and slab[a][3] > temp_high/2 : 
		outmost_al.append(slab[a])

print '''
(3) Collect the outmost layer with considering 
	(a) the number of broken bonds (less than ideal, maximun bonds, %5.i)
	(b) higher than middle of slab
	* since it is fcc metal system, the ideal nearest neighbor is '12'
	* the sub-surface is not considered, yet (will be updated)''' %max(N_bb)
# -----------------------------------------------------------------------------



top=[]; temp_t=[]
bridge=[]; temp_b=[]
hollow=[]; temp_h=[]
hollow4=[]; temp_h4=[]

temp_cell, temp_compound, expd_out_al = JH_lib.pst_cell_expansion(unitcell, outmost_al, dimension)
expd_out_al = JH_lib.pst_poscar_move(expd_out_al, float(unitcell[0][0])+float(unitcell[1][0]), float(unitcell[1][1]), 0)

for a in range(len(slab)):
	slab[a][0]=  slab[a][0].split('_')[0]

for a in outmost_al:
	for b in expd_out_al:
		if b != a and bonds_length * 0.9 < JH_lib.distance_atoms(a[1:4],b[1:4]) <  bonds_length * 1.1 : 	
			for c in expd_out_al:
				if c != a and bonds_length * 0.9 < JH_lib.distance_atoms(a[1:4],c[1:4]) < bonds_length * 1.1 and \
				   c != b and bonds_length * 0.9 < JH_lib.distance_atoms(b[1:4],c[1:4]) < bonds_length * 1.1 :
					for d in expd_out_al:
						if d != a and bonds_length * 0.9 < JH_lib.distance_atoms(a[1:4],d[1:4]) <  bonds_length * 1.1 and \
						   d != b and bonds_length * 0.9 < JH_lib.distance_atoms(b[1:4],d[1:4]) <  bonds_length * 1.1 and \
						   d != c and bonds_length * 0.9 < JH_lib.distance_atoms(c[1:4],d[1:4]) <  bonds_length * 1.1:
							# find 4 fold hollow site (4 bond)
							temp_h4.append(average([a[1:4],b[1:4],c[1:4],d[1:4]]))
					# find hollow site (3 bond)
					temp_h.append(average([a[1:4],b[1:4],c[1:4]]))
			# find bridge site (2 hand)
			if average([a[1:4],b[1:4]]) in temp_h4:	print 'in h4' ;pass
			else:	
				temp_b.append(average([a[1:4],b[1:4]]))
	# find top site (1 hand)	
	temp_t.append(a[1:4])

print '''
(4) Count the possible sites on the slab
	* the bridge site with overlapped to 4fH are neglected
	* Each bond-length of nearest atoms in slab from the chosen site are between %3.3f and %3.3f
''' %(bonds_length * 0.9 , bonds_length * 1.1 )

print '	(4-1) the top site    (T)	: have 1 bond  with surface (No. of maximum site: %3.i) ' %len(temp_t)
if len(temp_t)  != 0: 
	top=remove_overlap(temp_t, unitcell, adsorp)
	top = remove_overlap2_dist(adsorp, top)
print '\n	(4-2) the bridge site (B)	: have 2 bonds with surface (No. of maximum site: %3.i)' %len(temp_b)
if len(temp_b)  != 0: 
	bridge=remove_overlap(temp_b, unitcell, adsorp); 
	bridge = remove_overlap2_dist(adsorp, bridge)
print '\n	(4-3) the hollow site (H)	: have 3 bonds with surface (No. of maximum site: %3.i)' %len(temp_h)
if len(temp_h)  != 0: 
	hollow=remove_overlap(temp_h, unitcell, adsorp)
	hollow = remove_overlap2_dist(adsorp, hollow)
print '\n	(4-4) the 4-fold H  (4fH)	: have 4 bonds with surface (No. of maximum site: %3.i)' %len(temp_h4)
if len(temp_h4)  != 0: 
	hollow4=remove_overlap(temp_h4, unitcell, adsorp)
	hollow4 = remove_overlap2_dist(adsorp, hollow4)

# -----------------------------------------------------------------------------

print '''
(5) Generate the structure considering following factors
	(a) Possible position on surface (top, bridge, hollow, ...)
	(b) Interesting region between two circle having a radious of %4.3f and %4.3f, respectively
	(c) Default distance from the site to added adsorbate: %4.3f
	(d) for inclined/sloped part of surface, the normal vectors are calculated 
	(e) generated file names are site-name_integer.vasp (ex. top_1.vasp)
	\n\n\t\tWith reducing number of position, ''' \
%( 0, radious_from_adsroption, bond_length_finding(position) * 0.9 )

bond_length_finding(slab) #, bond_length_finding(outmost_al) 
# for a in range(len(slab)):
# 	for b in range(a):
# 		if JH_lib.distance_atoms(slab[a][1:4],slab[b][1:4]) <= bonds_length:
# 			print JH_lib.distance_atoms(slab[a][1:4],slab[b][1:4]), slab[a],slab[b]
#ex. generate_add_ads(site, radious_from_adsroption, orig_position, outmost_slab, unitcell, adsorp, name_of_site):
generate_add_ads(   top, radious_from_adsroption,      position,   outmost_al, unitcell, adsorp, 'top')
generate_add_ads(bridge, radious_from_adsroption,      position,   outmost_al, unitcell, adsorp, 'bridge')
generate_add_ads(hollow, radious_from_adsroption,      position,   outmost_al, unitcell, adsorp, 'hollow')
# # generate_add_ads(hollow,   [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'hollow')
# # generate_add_ads(hollow_4, [0,distance[0][0]+perterbation], slab, outmost_al,unitcell, adsorp, 'hollow4f')


# -----------------------------------------------------------------------------
print '''
(6) The generated files are stored in the folder of vasp
	(6-1) If there is folder named 'vasp' it will move as 'vasp_old'

# ----------------------------------- DONE ----------------------------------- #
#                                                                              #
#             YOU SHOULD CHECK THE GENERATED STRUCTURE ONE-BY-ONE              #
#             YOU SHOULD CHECK THE GENERATED STRUCTURE ONE-BY-ONE              #
#             YOU SHOULD CHECK THE GENERATED STRUCTURE ONE-BY-ONE              #
#                                                                              #
#                                                                              #
#  If you have any more question, Please send a mail to 'jkwnlee88@gmail.com'  #
#                                                                              #
# ----------------------------------- DONE ----------------------------------- #

''' 
os.system('if [[ `ls -l | grep vasp | awk \'{print $1}\'` ]]; then mv vasp vasp_old; mkdir vasp; fi')
os.system('mv *.vasp vasp/.')
