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
    1-1. post processing (pst)
        1-1-1. pst_poscar_move
        1-1-2. pst_poscar_ratation_atoms[not yet]
        1-1-3. pst_rot_compar_R[not yet]
        1-1-4. pst_cell_expansion
    1-2. molecule building
        1-2-1. mlc_ch4
        1-2-2. mlc_ch3
        1-2-3. mlc_water


2. add-ons (side-tools)
    2-1. poscar_head
         : Not arranged
    2-2. component_from_positions


3. judge_pri_conv
: For some structure, we prefer to have a lattice constant in conventional cell than primitive cell
  In order to do that, we use this function to get conventional cell lattice constant


4. w_poscar(filename, unitcell, compound, position, Selective):
: To write VASP crystal structure file type, this function need to include
(1) filename
(2) unitcell=[a, b, c] axis within xyz-coordinate unit of \AA
(3) compound=[[kind of compound],[number of atoms for each]]
(4) position=[[kind of compound, position of atom 1 in xyz], [kind of compound, position of atom 2 in xyz]...]
(5) Selective= True or False 
    This is required to apply Selective Dynamics in slab model or others


5. surface building: 
    5-1. scaling_making_slab(output_name, unitcell_print, initial_position,  \
                        test_x, test_y, diff_atom_position,       \
                        num_fix_lay, num_atom, seq_num, \
                        centering):

    5-2. build_slab(system, metal_kind, lat_con, index, output_name, \
                   vaccume,num_atom,num_fix_lay, symmetric):
            if system == 'fcc': 
                5-2-1. index = '100'
                5-2-2. index = '110'
                5-2-3. index = '111'
                5-2-4. index = '210'
                5-2-5. index = '211'
                5-2-6. index = '311'
                5-2-7. index = '331'
            if system == 'bcc': 'will be updated'

    5-3. add_adsorbate(slab_position, initial_position, mlc_position, output_name)
        # rotation will be updated

6. etc
    6-1 normal_vector(array_for_three_point)
        # Find a normal vector with three point
    6-2 cross_product(vec1,vec2)
        # Calculate cross product for two vector
    6-3 dot_product(vec1,vec2)
        # Calculate dot product for two vector
    6-4 shift_ads_w_normal_vec(inital_point, distance, vec)
        # Shift the initial point with following vector


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
    # poscar=[unitcell, compound, position]


## ------------------------------------------------------------------------------


def pst_poscar_move(position, mv_x, mv_y,mv_z):
    new_position = position
    for temp in range(len(position)):
        new_position[temp][1]=position[temp][1]-mv_x
        new_position[temp][2]=position[temp][2]-mv_y
        new_position[temp][3]=position[temp][3]-mv_z
    return new_position


## ------------------------------------------------------------------------------


def pst_cell_expansion(unitcell, compound, position, dimension):
    """
    1-1-4. pst_cell_expansion

    To make a supercell, expand the cell and add atoms properly
    output datas are new_unitcell, new_compound, new_position
    """
    new_position=[]; 
    for temp1 in range(len(position)):
        ## expanding to a
        for temp2 in range(dimension[0]): 
            ## expanding to b
            for temp3 in range(dimension[1]):
                ## expanding to c
                for temp4 in range(dimension[2]):
                    if len(position[temp1]) == 4: 
                        new_position.append(\
                            [position[temp1][0],\
                             position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
                             position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
                             position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , ])
                    if len(position[temp1]) == 7:
                        new_position.append(\
                            [position[temp1][0],\
                             position[temp1][1] + temp2 * float(unitcell[0][0]) + temp3 * float(unitcell[1][0]) + temp4 * float(unitcell[2][0]) , \
                             position[temp1][2] + temp2 * float(unitcell[0][1]) + temp3 * float(unitcell[1][1]) + temp4 * float(unitcell[2][1]) , \
                             position[temp1][3] + temp2 * float(unitcell[0][2]) + temp3 * float(unitcell[1][2]) + temp4 * float(unitcell[2][2]) , \
                             position[temp1][4], position[temp1][5], position[temp1][6]])

    new_unitcell = unitcell
    for temp1 in range(len(unitcell)):
        # each axis
        for temp2 in range(len(unitcell[temp1])):
            #each component
            new_unitcell[temp1][temp2] = float(unitcell[temp1][temp2]) * dimension[temp1]

    ## Count Component from positions
    new_compound = component_from_positions(new_position)

    return new_unitcell, new_compound, new_position


## ------------------------------------------------------------------------------


def mlc_ch4(initial_position_C=None, bond_leng=None):
    """
    1-2-1. mlc_ch4: building \ce{CH4} molecule 
    : As first varibalbe, define the position of C as like [0,0,0] 
      and 
      (Optional) insert the bonding length as second variable
                 [default is  1.09601 angstrom]

    : rotation will be [TODO:updated]

    """
    compound = [['C', 'H'], [1, 4]]
    if bond_leng == None:
        bond_leng = 1.09601
    if initial_position_C == None:
        initial_position_C = [0, 0, 0]

    position=[['C',                          0. ,                       0. ,              0. , 'T', 'T', 'T' ],\
              ['H',                          0. ,                       0. , -     bond_leng , 'T', 'T', 'T' ],\
              ['H',  8.**(0.5) * bond_leng / 3. ,                       0. ,   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. ,  bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. , -bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ] ]
    position = pst_poscar_move(position, - initial_position_C[0], - initial_position_C[1], - initial_position_C[2]) 

    return position, compound


## ------------------------------------------------------------------------------


def mlc_ch3(initial_position_C=None, bond_leng=None):
    """
    1-2-2. mlc_ch3: building \ce{CH3} molecule 
    : As first varibalbe, define the position of C as like [0,0,0] 
      and 
      (Optional) insert the bonding length as second variable
                 [default is  1.09601 angstrom]

    : rotation will be [TODO:updated]

    """
    compound = [['C', 'H'], [1, 3]]
    if bond_leng == None:
        bond_leng = 1.09601
    if initial_position_C == None:
        initial_position_C = [0, 0, 0]

    position=[['C',                          0. ,                       0. ,              0. , 'T', 'T', 'T' ],\
              ['H',  8.**(0.5) * bond_leng / 3. ,                       0. ,   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. ,  bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ],\
              ['H', bond_leng * -2.**(0.5) / 3. , -bond_leng * 2 / 6**(0.5),   bond_leng / 3 , 'T', 'T', 'T' ] ]
    position = pst_poscar_move(position, - initial_position_C[0], - initial_position_C[1], - initial_position_C[2]) 

    return position, compound


## ------------------------------------------------------------------------------


def poscar_head(num_atom,num_rlx):
    scale=( num_atom - 1 ) / 2
    num_fix_lay=( num_atom - num_rlx ) / 2 


## ------------------------------------------------------------------------------


def component_from_positions(insert_position):
    new_compound=[]
    componets = []
    num_componets = []
    i = 0; k= 0
    for temp1 in range(len(insert_position)):
        if temp1 == 0:
            componets.append(insert_position[temp1][0])
            num_componets.append(1)
        else:
            if insert_position[temp1][0] == temp:
                if componets[i] == insert_position[temp1][0]:
                    num_componets[i] = num_componets[i] + 1
            else:
                i = i + 1
                componets.append(insert_position[temp1][0])
                num_componets.append(1)

        temp = insert_position[temp1][0]


    new_compound = [componets , num_componets]
    return new_compound    

## ------------------------------------------------------------------------------


def judge_pri_conv(filename, system):
    ## Function check the condition of file wheather it is primitive or not
    ##
    ## Condition1: number atom == 1
    ## Condition2: nit
    unitcell, compound, position=r_cryst_vasp('CONTCAR')
    if len(compound[1]) == 1:
        # print "this cell include only one kind of atom"
        if compound[1][0] == '1':
            # print "this cell include only one of atom"
            True
    else:
        print "this cell include only more than one kind of atoms", compound[1]
        pass
    if system == 'fcc':
        material_name=compound[0][0]
        latcon=unitcell[0][0]
        if float(latcon) == 0.:
            # print "the unit cell is not conventional"
            latcon=float(unitcell[0][1]) * 2
        else:
            pass
    return latcon


## ------------------------------------------------------------------------------


def w_poscar(compound, position, filename = None, unitcell = None, Selective = None):
    """
    Write a poscar file type
    Officially, it needs 2 parameters. "compound and position"

    [default] filename is 'py_POSCAR.vasp'
    [default] unitcell is a cubic with a lattice, 1.1 times bigger than maximum in "position"
    [default] Selective == None, if you want to turn on, 
              Please type True with proper "position" including Selective information inside
              ( ex. [[Cu, 0, 0, 0, F, F, F]]  )
    """
    if filename == None:
        filename = 'py_POSCAR.vasp'

    if unitcell == None:
        tempx = 10.
        tempy = 10.
        tempz = 10.
        for temp in position:
            tempx = max(abs(float(temp[1])), tempx)
            tempy = max(abs(float(temp[2])), tempy)
            tempz = max(abs(float(temp[3])), tempz)
        unitcell = [[tempx*1.1, 0., 0.], \
                    [ 0.,tempy*1.1, 0.], \
                    [ 0., 0.,tempz*1.1]]

    if len(position[1]) == 7:
        Selective = True


    f=open(filename, 'w')
    f.write('test\n %19.16f\n' % 1.00000000000000)
    for temp in range(len(unitcell)):
        xyz=unitcell[temp]
        if len(xyz) == 3:
            f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print " WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS"
            break

    for temp in range(len(compound)):
        for temp2 in range(len(compound[temp])):
            f.write(' %4.3s' %compound[temp][temp2])
        f.write('\n')
    if Selective == None:
        pass
    else:
        f.write('Selective Dynamics\n')
    f.write('C\n')

    for temp in range(len(position)):
        xyz=position[temp][1:4]
        if len(xyz) == 3:
            if Selective:
                dynamic=position[temp][4:]
                f.write(' %22.16f%22.16f%22.16f %5.1s %5.1s %5.1s \n' \
                    % (float(xyz[0]), float(xyz[1]), float(xyz[2]), dynamic[0],dynamic[1], dynamic[2])) 
            else:
                f.write(' %22.16f%22.16f%22.16f \n' % (float(xyz[0]), float(xyz[1]), float(xyz[2]))) 
        else: 
            print " WARNING the xyz is not proper IT SHOULD CONTAIN ONLY 3 NUMBERS"
            break
## TEST for w_poscar
# unitcell, compound, positions=r_cryst_vasp('CONTCAR')
# positions[0].append('T')
# positions[0].append('T')
# positions[0].append('T')
# w_poscar('test.vasp', unitcell, compound, positions, False)


## ------------------------------------------------------------------------------


def scaling_making_slab(output_name, unitcell_print, initial_position,  \
                        test_x, test_y, diff_atom_position,       \
                        num_fix_lay, num_atom, seq_num, \
                        centering):
    ## store positions into positions
    position_temp=initial_position[0][1:4]
    position_after=position_temp
    position_print=[]
    for temp in range(num_atom):
        position_print.append(['Cu',position_temp[0], position_temp[1], position_temp[2]])
        ## Define x (a) coefficient!
        position_after[0] = position_temp[0] + diff_atom_position[0]    
        if position_after[0] >= test_x: position_after[0]= position_after[0] - test_x
        ## Define y (b) coefficient!
        position_after[1] = position_temp[1] + diff_atom_position[1]    
        if position_after[1] >= test_y: position_after[1]= position_after[1] - test_y
        ## Define z (c) coefficient!
        position_after[2] = position_temp[2] + diff_atom_position[2]
        position_temp=position_after

    if centering:
        z_collect=[]
        temp_num_rlx=1 ; temp_num_unrlx=1
        for temp in range(len(position_print)):
            z_collect.append(position_print[temp][3])
        averaged=reduce(lambda x, y: x + y, z_collect) / len(z_collect)
        for temp in range(len(position_print)):
            position_print[temp][3]=position_print[temp][3]-averaged
            if ( num_atom - num_fix_lay ) / 2.  >=  temp_num_rlx :
                position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')
                temp_num_rlx=temp_num_rlx + 1
            else:
                if temp_num_unrlx > num_fix_lay:
                    position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')
                else:
                    position_print[temp].append('F');position_print[temp].append('F');position_print[temp].append('F')
                    temp_num_unrlx = temp_num_unrlx + 1
    else:
        for temp in range(num_atom):
            if temp <=  num_fix_lay :
                position_print[temp].append('F');position_print[temp].append('F');position_print[temp].append('F')
            else:
                position_print[temp].append('T');position_print[temp].append('T');position_print[temp].append('T')

    compound_print=[[position_print[0][0]],[len(position_print)]]
    w_poscar(compound_print, position_print, output_name, unitcell_print, True)


## ------------------------------------------------------------------------------


def build_slab(system, metal_kind, lat_con, index, output_name, vaccume,num_atom,num_fix_lay, symmetric):
    if system == 'fcc':
        ## CHOOSE face-centered crystal
        position=[[metal_kind, 0.0, 0.0, 0.0, 'F', 'F', 'F']]
        if index == '100':
            unitcell_print, diff_atom_position = poscar_fcc_100_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=2
        elif index == '110':
            unitcell_print, diff_atom_position = poscar_fcc_110_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=2
        elif index == '111':
            unitcell_print, diff_atom_position = poscar_fcc_111_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=3
        elif index == '210':
            unitcell_print, diff_atom_position = poscar_fcc_210_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=10
        elif index == '211':
            unitcell_print, diff_atom_position = poscar_fcc_211_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=6
        elif index == '311':
            unitcell_print, diff_atom_position = poscar_fcc_311_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=11
        elif index == '331':
            unitcell_print, diff_atom_position = poscar_fcc_331_slab_primi(metal_kind, lat_con, vaccume,num_atom)
            seq_num=19

    test_x  = unitcell_print[0][0] + unitcell_print[1][0] + unitcell_print[2][0]
    test_y  = unitcell_print[0][1] + unitcell_print[1][1] + unitcell_print[2][1]

    scaling_making_slab(output_name, unitcell_print, position, \
        test_x ,test_y ,diff_atom_position,\
        num_fix_lay,  num_atom, seq_num, symmetric)


## ------------------------------------------------------------------------------


def poscar_fcc_100_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con / 2.**(0.5) ;vec_ay = 0.0 ;vec_az = 0.0
    vec_bx = 0.0         ;vec_by = lat_con / 2.**(0.5) ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 2 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]
    ## Position of next atom 
    len_x = vec_ax / 2.  ; len_y= vec_by / 2.  ;len_z = lat_con / 2.
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position     


## ------------------------------------------------------------------------------


def poscar_fcc_110_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con         ;vec_ay = 0.0         ;vec_az = 0.0
    vec_bx = 0.0             ;vec_by = lat_con / 2.**(0.5)  ;vec_bz = 0.0
    vec_cx = 0.0             ;vec_cy = 0.0
    vec_cz = lat_con / 2.**(0.5) / 2 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]
    ## Position of next atom 
    len_x = vec_ax  / 2.    ;len_y = vec_by  / 2.  ;len_z = lat_con / 2.**(0.5) / 2
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_111_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con / 2.**(0.5)  ;vec_ay = 0.0     ;vec_az = 0.0
    vec_bx = lat_con / 2.**(0.5) / 2.
    vec_by = lat_con / 2.**(0.5) / 2. * 3. **(0.5) ;vec_bz = 0.0
    vec_cx = 0.0                  ;vec_cy = 0.0
    vec_cz = lat_con * 3. **(0.5) / 3 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom 
    len_x = ( vec_bx + vec_ax) / 3. ; len_y = vec_by / 3. ;len_z = lat_con * 3. **(0.5) / 3
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_211_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 2.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = 0.0         ;vec_by = lat_con * 3.**(0.5)     ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con * 6. **(0.5) / 2 / 6 * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom 
    len_x = vec_ax / 2.  ;len_y = vec_by * 2. / 3.         ;len_z = lat_con * 6. **(0.5) / 2. / 6.
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_311_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 6.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = lat_con * 6.**(0.5) / 2. * 5 / 6. 
    vec_by = lat_con * 6.**(0.5) / 2  * 11.**(0.5) / 6     ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 11. **(0.5) * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * 0.2727  ;len_y = vec_by * 0.2727  ;len_z = lat_con / 11. **(0.5)
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_210_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 6.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = vec_ax * 2. / 3.  ;vec_by = vec_ax * 5.**(0.5) / 3  ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con * 5.**(0.5) / 10. * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * 0.7  ;len_y = vec_by * 0.7  ;len_z = lat_con * 5. **(0.5) / 10
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def poscar_fcc_331_slab_primi(metal_kind, lat_con, vaccume,num_atom):
    ## Build Periodic Boundary Condition
    vec_ax = lat_con * 10.**(0.5) / 2. ;vec_ay = 0.0        ;vec_az = 0.0
    vec_bx = lat_con * 10.**(0.5) / 2. * 0.9
    vec_by = lat_con * 10.**(0.5) / 2. * (0.19)**(0.5)      ;vec_bz = 0.0
    vec_cx = 0.0         ;vec_cy = 0.0
    vec_cz = lat_con / 19.**(0.5)  * ( num_atom - 1 ) + vaccume
    unitcell_print=[[vec_ax, vec_ay, vec_az], [vec_bx, vec_by, vec_bz], [vec_cx, vec_cy, vec_cz]]

    ## Position of next atom  (CONTROVERSY for 0.2727)
    len_x = (vec_ax + vec_bx ) * (1 - 10.**(0.5) / 10 )
    len_y = vec_by * (1 - 10.**(0.5) / 10 )  ;len_z = lat_con / 19. **(0.5)
    diff_atom_position = [len_x, len_y, len_z]
    return unitcell_print, diff_atom_position


## ------------------------------------------------------------------------------


def add_adsorbate(slab_position, slab_unitcell, adsorbate_position, output_name):
    new_position = slab_position + adsorbate_position
    new_compound = component_from_positions(new_position)
    # print new_position, new_compound
    w_poscar(new_compound, new_position, output_name, slab_unitcell)


## ------------------------------------------------------------------------------


def print_error(True):
    print """
 _____ ____  ____   ___  ____
| ____|  _ \|  _ \ / _ \|  _ \\
|  _| | |_) | |_) | | | | |_) |
| |___|  _ <|  _ <| |_| |  _ <
|_____|_| \_\_| \_\\\___/|_| \_\\
"""

## ------------------------------------------------------------------------------


def normal_vector(array_for_three_point):
    """
    Calculate the normal vector with a inserted three points
    """

    if len(array_for_three_point) != 3:
        return print_error(True)
    point0 = array_for_three_point[0]; point1 = array_for_three_point[1];point2 = array_for_three_point[2]
    vec1=[]
    vec2=[]
    center_of_points=[]
    for a in range(3):
        vec1.append(point0[a] - point1[a])
        vec2.append(point0[a] - point2[a])
        center_of_points.append( (point0[a]+ point1[a] + point2[a]) / 3. )
    print vec1
    print vec2
    if cross_product(vec1, vec2)[2] < 0:
        return cross_product(vec2, vec1), center_of_points
    else:
        return cross_product(vec1, vec2), center_of_points
    # if cross_vec12


## ------------------------------------------------------------------------------


def cross_product(vec1, vec2):
    """
    Calculate the corss product from vec1 to vec2
    """
    if len(vec1) != 3 or len(vec2) != 3:
        return print_error(True)
    else:
        return [vec1[1] * vec2[2]  - vec1[2] * vec2[1], 
                vec1[2] * vec2[0]  - vec1[0] * vec2[2],
                vec1[0] * vec2[1]  - vec1[1] * vec2[0]]


## ------------------------------------------------------------------------------


def shift_ads_w_normal_vec(inital_point,distance, vec):
    """
    Shift the initial point with following vector
    """
    if len(inital_point) !=3 or len(vec) !=3:
        return print_error(True)
        
    len_vec = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    temp_vec = []; next_position=[]
    for temp in range(3):
        temp_vec.append(vec[temp]/len_vec * distance)
        next_position.append(inital_point[temp] + temp_vec[temp])
    print(' %22.16f%22.16f%22.16f T  T  T\n' % (float(next_position[0]), float(next_position[1]), float(next_position[2]))) 

    return temp_vec,next_position


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
# unitcell, compound, position=r_cryst_vasp('POSCAR')
scale=5
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
    if a[0] == 'Zr':
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
        system.append(['Zr', x, y, z])

for a in system:
    for b in range(len(grid)):
        tarbulate[b] = tarbulate[b] + counts(a[1], a[2],a[3], grid[b][0], grid[b][1], grid[b][2], scale)

i = 0
print xwmax, ywmax, zwmax
for a in tarbulate:
    print "%15.10f" %(a) ,
    if i == 4:
        i = 0
        print ""
    else:
        i = i + 1









