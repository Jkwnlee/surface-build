#!/bin/python

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


4. write structure
    4.1 w_poscar(filename, unitcell, compound, position, Selective):
    : To write VASP crystal structure file type, this function need to include
    (1) filename
    (2) unitcell=[a, b, c] axis within xyz-coordinate unit of \AA
    (3) compound=[[kind of compound],[number of atoms for each]]
    (4) position=[[kind of compound, position of atom 1 in xyz], [kind of compound, position of atom 2 in xyz]...]
    (5) Selective= True or False 
        This is required to apply Selective Dynamics in slab model or others

    4.2 w_lammps_str(filename, unitcell, compound, position, Selective):
    : To write LAMMPS crystal structure file type, 
    4.2.1 mass_of_element(compound)


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
    4.1 Write a poscar file type
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


## ------------------------------------------------------------------------------


def w_lammps_str(compound, position, filename = None, unitcell = None, sequence_of_compound = None):
    """
    4.2 w_lammps_str(filename, unitcell, compound, position, Selective):
    : To write LAMMPS crystal structure file type, 
    Officially, it needs 2 parameters. "compound and position"

    [default] filename is 'py_POSCAR.vasp'
    [default] unitcell is a cubic with a lattice, 1.1 times bigger than maximum in "position"
    [default] Selective == None, if you want to turn on, 
              Please type True with proper "position" including Selective information inside
              ( ex. [[Cu, 0, 0, 0, F, F, F]]  )
    """
    if filename == None:
        filename = 'py_structure.data'

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

    if sequence_of_compound == None:
        sequence_of_compound = compound[0]  ## TODO later
    n_atom=0
    for a in compound[1]:
        n_atom=n_atom + int(a)

    f=open(filename, 'w')
    f.write('# LAMMPS data file created from JH\'s python\n')
    f.write('%5.0f atoms\n' %(n_atom ) )
    f.write('%3.0f atom types\n' %(len(compound[0]) ) )
    f.write('%10.6f %8.6f  xlo xhi\n' %( 0, float(unitcell[0][0]) ) )
    f.write('%10.6f %8.6f  ylo yhi\n' %( 0, float(unitcell[1][1]) ) )
    f.write('%10.6f %8.6f  zlo zhi\n' %( 0, float(unitcell[2][2]) ) )
    f.write('%10.6f %8.6f %8.6f  xy xz yz\n\nMasses\n\n' %( 0,0,0  ) )

    compound = mass_of_element(compound)
    i = 1
    for a in sequence_of_compound[0]:
        for b in range(len(compound[0])):
            if a == compound[0][b]:
                f.write('%1.0f %8.6f\n' %(i, float(compound[2][b])))
                i = i + 1
    f.write('\nAtoms\n\n')

    i = 1
    for a in range(len(sequence_of_compound[0])):    
        for temp in position:
            if sequence_of_compound[0][a] == temp[0]:
                f.write('%1.0f %1.0f %2.2s %8.6f %8.6f %8.6f\n' \
                         %(i, a+1, sequence_of_compound[1][a], float(temp[1]), float(temp[2]), float(temp[3]) ) )
                i = i + 1


## ------------------------------------------------------------------------------


def mass_of_element(compound):
    lib=[['1.0079','H'], ['4.0026','He'], ['6.941','Li'], ['9.0122','Be'],\
        ['10.811','B'], ['12.0107','C'], ['14.0067','N'], ['15.9994','O'],\
        ['18.9984','F'], ['20.1797','Ne'], ['22.9897','Na'], ['24.305','Mg'],\
        ['26.9815','Al'], ['28.0855','Si'], ['30.9738','P'], ['32.065','S'],\
        ['35.453','Cl'], ['39.0983','K'], ['39.948','Ar'], ['40.078','Ca'],\
        ['44.9559','Sc'], ['47.867','Ti'], ['50.9415','V'], ['51.9961','Cr'],\
        ['54.938','Mn'], ['55.845','Fe'], ['58.6934','Ni'], ['58.9332','Co'],\
        ['63.546','Cu'], ['65.39','Zn'], ['69.723','Ga'], ['72.64','Ge'],\
        ['74.9216','As'], ['78.96','Se'], ['79.904','Br'], ['83.8','Kr'],\
        ['85.4678','Rb'], ['87.62','Sr'], ['88.9059','Y'], ['91.224','Zr'],\
        ['92.9064','Nb'], ['95.94','Mo'], ['98','Tc'], ['101.07','Ru'],\
        ['102.9055','Rh'], ['106.42','Pd'], ['107.8682','Ag'], ['112.411','Cd'],\
        ['114.818','In'], ['118.71','Sn'], ['121.76','Sb'], ['126.9045','I'],\
        ['127.6','Te'], ['131.293','Xe'], ['132.9055','Cs'], ['137.327','Ba'],\
        ['138.9055','La'], ['140.116','Ce'], ['140.9077','Pr'], ['144.24','Nd'],\
        ['145','Pm'], ['150.36','Sm'], ['151.964','Eu'], ['157.25','Gd'],\
        ['158.9253','Tb'], ['162.5','Dy'], ['164.9303','Ho'], ['167.259','Er'],\
        ['168.9342','Tm'], ['173.04','Yb'], ['174.967','Lu'], ['178.49','Hf'],\
        ['180.9479','Ta'], ['183.84','W'], ['186.207','Re'], ['190.23','Os'],\
        ['192.217','Ir'], ['195.078','Pt'], ['196.9665','Au'], ['200.59','Hg'],\
        ['204.3833','Tl'], ['207.2','Pb'], ['208.9804','Bi'], ['209','Po'],\
        ['210','At'], ['222','Rn'], ['223','Fr'], ['226','Ra'], ['227','Ac'],\
        ['231.0359','Pa'], ['232.0381','Th'], ['237','Np'], ['238.0289','U'],\
        ['243','Am'], ['244','Pu'], ['247','Cm'], ['247','Bk'], ['251','Cf'],\
        ['252','Es'], ['257','Fm'], ['258','Md'], ['259','No'], ['261','Rf'],\
        ['262','Lr'], ['262','Db'], ['264','Bh'], ['266','Sg'], ['268','Mt'],\
        ['272','Rg'], ['277','Hs'], ['Darmstadtium','110'], ['Ununbium','112'],\
        ['Ununtrium','113'], ['Ununquadium','114'], ['Ununpentium','115'], ['Ununhexium','116'],\
        ['Ununseptium','117'], ['Ununoctium','118']]
    temp=[]
    for a in compound[0]:
        for b in lib:
            if a == b[1]:
                temp.append(b[0])

    compound.append(temp)
    return(compound)


## ------------------------------------------------------------------------------

# target = str(sys.argv[1])
target ='c00016.vasp'
sequence_of_compound = [['Zr', 'O', 'Y'],['+4','-2','+3']]
unitcell, compound, position = r_cryst_vasp(target)
w_lammps_str(compound, position, "structure.data", unitcell, sequence_of_compound)





# ###########################################################################
# ###########################################################################
# ####   ____   ___    _   _  ___ _____   _____ ___  _   _  ____ _   _   ####
# ####  |  _ \ / _ \  | \ | |/ _ \_   _| |_   _/ _ \| | | |/ ___| | | |  ####
# ####  | | | | | | | |  \| | | | || |     | || | | | | | | |   | |_| |  ####
# ####  | |_| | |_| | | |\  | |_| || |     | || |_| | |_| | |___|  _  |  ####
# ####  |____/ \___/  |_| \_|\___/ |_|     |_| \___/ \___/ \____|_| |_|  ####
# ####                                                                   ####
# ###########################################################################
# ###########################################################################
