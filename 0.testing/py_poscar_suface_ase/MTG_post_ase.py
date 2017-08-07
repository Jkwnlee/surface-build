#!/bin/python

## by Ji-Hwan Lee (jkwnleee88@gmail.com)
## 23th Feb 2017
## purpose

import math

### Number of each atoms
def number_of_atom_kinds():
  with open('POSCAR', 'r') as cont:
    lines = cont.readlines()
    tmp = lines[5]
    kind_of_atoms = tmp.split()
    return kind_of_atoms

def RepresentsFlt(s):
  try: 
      float(s)
      return True
  except ValueError:
      return False

def check_direct(let):
  if let[0] == 'D' or let[0] == 'd':return 'Direct'
  elif let[0] == 'S': return 'Selective'
  elif let[0] == 'C' or 'c': return 'Cartesian'
  # else: return 'pass'


def read_poscar(filename):
  head=[] 
  scale=[]
  vector=[]
  kind_of_atoms=[]
  num_of_atoms=[]
  system=[]
  atom_position=[]
  with open(filename, 'r') as posc:
    i=0
    for line in posc:
      if i == 0: head=line.split()
      elif i == 1: scale=line.split()
      elif i == 2 or i == 3 or i == 4: vector.append(line.split()[:] )
      elif i == 5: 
        if RepresentsFlt(line.split()[0]): 
          num_of_atoms=line.split()[0]
          kind_of_atoms=head
        else: kind_of_atoms=line.split()[0]
      elif i==6:
        if num_of_atoms == []: num_of_atoms=line.split()[0]
        else: system = check_direct(line.split()[0])
      elif i== 7:
        if system == []: system = check_direct(line.split()[0])
        else: pass
      elif i == 8:
        if system == 'Selective': system = check_direct(line.split()[0])
        else: pass
      else: atom_position.append(line.split())
      i= i+1

  # print head
  # print scale
  # print vector 
  # print kind_of_atoms
  # print num_of_atoms
  # print system
  # print atom_position
  return {'head':head,
    'scale':scale,
    'vector':vector,
    'kind_of_atoms':kind_of_atoms,
    'num_of_atoms':num_of_atoms,
    'system':system,
    'atom_position':atom_position
    }


def find_bigger(reference, target):
  if target == reference or target > reference:
    bigger = target
  elif target < reference:
    bigger = reference  

  return bigger

def post_ase_vasp_rewrite(filename,target_filename,rearrange,  vacuum_distance,  select_dynamic):
  head=[] 
  scale=[]
  vector=[]
  kind_of_atoms=[]
  num_of_atoms=[]
  system=[]
  atom_position=[]
  with open(filename, 'r') as posc:
    i=1
    for line in posc:
      if i == 1: head=line.split()[0]
      elif i == 2: scale=line.split()[0]
      elif i == 3 or i == 4 or i == 5: vector.append(line.split()[:] )
      elif i == 6: 
        if RepresentsFlt(line.split()[0]): 
          num_of_atoms=line.split()[0]
          kind_of_atoms=head
        else: kind_of_atoms=line.split()[0]
      elif i==7:
        if num_of_atoms == []: num_of_atoms=line.split()[0]
        else: system = check_direct(line.split()[0])
      elif i==8:
        if system == []: system = check_direct(line.split()[0])
        else: atom_position.append(line.split())
      elif i == 9:
        if system == 'Selective': system = check_direct(line.split()[0])
        else: atom_position.append(line.split())

      else: atom_position.append(line.split())
      i= i+1



  with open(target_filename, 'w') as posc:
    #######################################
    ## START TO STORE THE ATOMS POSITION ##
    #######################################
    data_m=[]
    ## Start to rearrange atom position [ing]
    i=0; highest_z=0.
    if rearrange =='True': 
      ## Collect z value and calculate average z: temp_val_avez
      temp_val_avez = float(vector[2][2]) /2
      
      ## Replace the z position 
      for temp in range(len(atom_position)):
        x=[]
        x=atom_position[temp]
        vec_x=x[0];vec_y=x[1];vec_z=x[2]
        temp_val_z = float(vec_z) - temp_val_avez
        x[2]=temp_val_z 
        if abs(temp_val_z) < select_dynamic[2] / 2:
          if select_dynamic[0] == 'True' :
            if abs(temp_val_z) > select_dynamic[1] / 2:
              data="%24.15f"%float(x[0]) + "%24.15f"%float(x[1]) + "%24.15f"%float(x[2]) + " T T T"
              data_m.append(data)
              i = i+1; highest_z = find_bigger(highest_z, abs(float(x[2])) )
            else:
              data="%24.15f"%float(x[0]) + "%24.15f"%float(x[1]) + "%24.15f"%float(x[2]) + " F F F"
              data_m.append(data)
              i = i+1; highest_z = find_bigger(highest_z, abs(float(x[2])) )
          else:
            data="%24.15f"%float(x[0]) + "%24.15f"%float(x[1]) + "%24.15f"%float(x[2])
            data_m.append(data)
            i = i+1; highest_z = find_bigger(highest_z, abs(float(x[2])) )
          # posc.write(data+'\n')
          data_m.append('\n')
        else:
          pass

    ## If the rearrangement is not required [DONE]
    else:
      for temp in range(len(atom_position)):
        x=[]
        x=atom_position[temp]
        for temp2 in range(len(x)):
          # posc.write("%24.15f"%float(x[temp2]))
          data_m.append("%24.15f"%float(x[temp2]))
          i = i+1; highest_z = find_bigger(highest_z, abs(float(x[2])) )
        # posc.write('\n')
        data_m.append('\n')


    #######################################
    ## START TO STORE THE HEAD OF SYSTEM ##
    #######################################

    data_mx=[] # matrix to store data
    ## Write a head, scale of POSCAR
    data_mx.append(head+'\n'); data_mx.append(scale+'\n')
      ## Write a unitcell vector, 
    for temp in range(len(vector)):
      x=[]
      x=vector[temp]   
      if temp == 2:
        if i == num_of_atoms:
          pass
        else:
          num_of_atoms = str(i)
          data="%24.18f"%float(x[0]) + "%24.18f"%float(x[1]) + "%24.18f"%float(highest_z * 2 +  vacuum_distance )
      else:
        data="%24.18f"%float(x[0]) + "%24.18f"%float(x[1]) + "%24.18f"%float(x[2])
      data_mx.append(data+'\n')
      # posc.write(data+'\n')
    if i == num_of_atoms:
      pass
    else:
      num_of_atoms = str(i)

    if select_dynamic[0] == 'True':
      ## Write a kind of atom, number of atom, selective dynamics (choose)
      # posc.write(kind_of_atoms+'\n'+num_of_atoms+'\nSelective dynamics\n'+system+'\n')
      data_mx.append(kind_of_atoms+'\n'+num_of_atoms+'\nSelective dynamics\n'+system+'\n')
    else:
      ## Write a kind of atom, number of atom
      # posc.write(kind_of_atoms+'\n'+num_of_atoms+'\n'+system+'\n')
      data_mx.append(kind_of_atoms+'\n'+num_of_atoms+'\n'+system+'\n')


    #######################################
    ## START TO WRITE THE DATA OF SYSTEM ##
    #######################################
    for length_data in range(len(data_mx)): posc.write(data_mx[length_data])

    for length_data in range(len(data_m)): posc.write(data_m[length_data])


  return True
