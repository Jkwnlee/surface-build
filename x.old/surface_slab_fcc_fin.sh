# INPUT: FCC optimized CONTCAR files (Output of VASP file)

######################################################################
#                             HOW to USE                             # 
#
# function_name $num_for_vac $num_for_atomes_in_symmetric_slab	$num_of_allowing_relaxation $asym (if you want to make asym type 1)
# function name list
# poscar_fcc_111_slab_primi
# poscar_fcc_110_slab_primi
# poscar_fcc_100_slab_primi
# poscar_fcc_210_slab_primi
# poscar_fcc_211_slab_primi
# poscar_fcc_311_slab_primi
# poscar_fcc_331_slab_primi
# poscar_fcc_557_slab_primi
#
# example: poscar_fcc_557_slab_primi 16 31 16
# Make a of 557 oriented FCC slab having 31 atomic layers with 16 Å of Vacuum and allowing relaxation for outmost 16 atomic layer 
#
######################################################################

######################################################################
#                      Recently updated contents                     #
# Find unit cell and repeating layer vacuum / slab / vaccum (not perfect just unit cell! )
# A. Find unit vector of unit cell *(Using vesta > Lattice plane)
# B. Define vector a/|a| and b/|b| ex) a/|a|=(0.5 0.5 0) b/|b|=(-0.5 3 2.5)
# C. Based on dot product, generate the a and b, ex) bx= cos(theta) b/|b| = 5/2/sqrt(31)
# D. Define vector  between original and atom in second layer (based on your guessing) ex)for FCC(557), z=(0 0.5 0.5) 
# E. The number of "len_y" is defined by  (total number of unit cell repeating, each step:-1)  * vec_by / total number of atom to repeat
#
# 20150319 Update: 557 surface
# 20150331 Update: asym
# 20150414 Update: primitive / conventional ->> latice constant
# 20150414 Update: considering tau_hkl as input for strained system  (tau_hkl comes from gibbs2 datas or E-z-strain datas)
# 20150414 Update: strianed commend for 111 direction
######################################################################

function poscar_head(){
	num_atom=$1; 	num_rlx=$2
	scale=$( echo "scale=30 ; ( $num_atom - 1 ) / 2  " | bc )
	num_fix_lay=$(echo "scale=30 ; ( $num_atom - $num_rlx ) / 2 " | bc |awk '{printf "%2.0f", $1}')
}

function poscar_lattic_from_contcar_prim(){
	material_name=$(cat ./CONTCAR |head -6 |tail -1| awk '{print $1}')
	scale_num=$(cat ./CONTCAR |head -2 | tail -1)
	lat_con_a=$(cat ./CONTCAR |head -3 | tail -1 | awk '{printf "%20.9f" , $1}')
	lat_con=$(echo "scale=30 ; $lat_con_a * $scale_num " | bc )
	if [[ `echo $lat_con_a |cut -c 1-3` == "0.0" ]]; then
# echo eheradingyo previous $lat_con_a
		lat_con_a=$(cat ./CONTCAR |head -`expr 3`| tail -1 | awk '{printf "%20.9f" , $2}')
		lat_con=$(echo "scale=30 ; 2 * $lat_con_a * $scale_num " | bc )
	fi
# echo 1 $material_name 2 $scale_num 3 $lat_con_a 4 $lat_con
}

function poscar_slab_vector_selec(){
	ax=$1;	bx=$2; 	by=$3;	cz=$4
	echo $ax"         0.0000000000000         0.0000000000000" >> ./POSCAR
	echo $bx"         "$by"         0.0000000000000" >> ./POSCAR
	echo "0.0000000000000         0.0000000000000         "$cz >> ./POSCAR
	echo $material_name >> ./POSCAR; 	    echo $num_atom >> ./POSCAR
	echo Selective Dynamics >> ./POSCAR;	echo Cartesian >> ./POSCAR
	echo "0.0000000000000         0.0000000000000         0.0000000000000  F F F" >> ./POSCAR  ## A layer
}


function scaling_making_slab(){
##	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	test_x=$1; 	test_y=$2; 	len_x=$3; 	len_y=$4; 	len_z=$5; num_fix_lay=$6 ; seq_num=$7 

	# cb=0; 	ca=0;	cz=0;
	cb2=0; 	ca2=0;	cz2=0;
	for a in `seq 1 1 $scale`
		do
		left=$( echo  $a%$seq_num |bc )
##Define x (a) coefficient!
		ca1=$( echo "scale=30 ; $ca2 + $len_x " |bc )
		if [ $(bc <<< "$ca1 < $test_x") -eq 1  ]; then
			ca2=$( echo "$ca1" |cut -c 1-15)             ## WHEN ca1  < $test
			else
			ca2=$( echo "scale=30 ; $ca1 - $test_x " |bc ) ## WHEN ca1  < $test
		fi
		ca=$( echo "$ca2" |cut -c 1-15)             ## WHEN ca1  < $test
		caa=$( echo "scale=30 ; $test_x - $ca  " |bc |cut -c 1-15)
##Define y (b) coefficient!
		cb1=$( echo "scale=30 ; $cb2 + $len_y" |bc )
		if [ $(bc <<< "$cb1 < $test_y") -eq 1 ]; then
			cb2=$( echo "$cb1" )             ## WHEN cb1  < $test
			else
			cb2=$( echo "scale=30 ; $cb1 - $test_y " |bc )
		fi
		cb=$( echo "$cb2" |cut -c 1-15)
		cbb=$( echo "scale=30 ; $test_y - $cb  " |bc |cut -c 1-15)
		if [ $left -eq 0 ]; then
			# echo $left
			ca=0.0000000000000
			cb=0.0000000000000
			cbb=$( echo "scale=30 ; $test_y - $cb  " |bc |cut -c 1-15)
			caa=$( echo "scale=30 ; $test_x - $ca  " |bc |cut -c 1-15)
		fi
##Define y (c) coefficient!
		cz1=$( echo "scale=30 ; $cz2 + $len_z " |bc |cut -c 1-15)
		cz=$( echo "scale=30 ; $cz1"  | bc | cut -c 1-15 )
		czz=$( echo "scale=30 ; $vec_cz - $cz" |bc |cut -c 1-15)
##Print the result!
		if [ $a -le  $num_fix_lay ] ;then
			echo $ca"         "$cb"         "$cz" F F F " >> ./POSCAR
			echo $caa"         "$cbb"         "$czz" F F F" >> ./POSCAR
			else
			echo $ca"         "$cb"         "$cz" T T T" >> ./POSCAR
			echo $caa"         "$cbb"         "$czz" T T T" >> ./POSCAR
		fi
	done
}



function scaling_making_asym_slab(){
##	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	# echo "# PLEAE check this poscar before using, the vacuum and something is not optimized" >> ./POSCAR
	test_x=$1; 	test_y=$2; 	len_x=$3; 	len_y=$4; 	len_z=$5; num_fix_lay=$6 ; seq_num=$7

	cb=0; 	ca=0;	cz=0
	for a in `seq 1 1 $scale`
		do
		left=$( echo  $a%$seq_num |bc )
##Define x (a) coefficient!
		ca1=$( echo "scale=30 ; $ca + $len_x " |bc )
		if [ $(bc <<< "$ca1 < $test_x") -eq 1  ]; then
			ca2=$( echo "$ca1" |cut -c 1-15)             ## WHEN ca1  < $test
			else
			ca2=$( echo "scale=30 ; $ca1 - $test_x " |bc ) ## WHEN ca1  < $test
		fi
		ca=$( echo "$ca2" |cut -c 1-15)             ## WHEN ca1  < $test
		caa=$( echo "scale=30 ; $test_x - $ca  " |bc |cut -c 1-15)
##Define y (b) coefficient!
		cb1=$( echo "scale=30 ; $cb + $len_y" |bc )
		if [ $(bc <<< "$cb1 < $test_y") -eq 1 ]; then
			cb2=$( echo "$cb1" )             ## WHEN cb1  < $test
			else
			cb2=$( echo "scale=30 ; $cb1 - $test_y " |bc )
		fi
		cb=$( echo "$cb2" |cut -c 1-15)
		cbb=$( echo "scale=30 ; $test_y - $cb  " |bc |cut -c 1-15)
		if [ $left -eq 0 ]; then
			# echo $left
			ca=0.0000000000000
			cb=0.0000000000000
			cbb=$( echo "scale=30 ; $test_y - $cb  " |bc |cut -c 1-15)
			caa=$( echo "scale=30 ; $test_x - $ca  " |bc |cut -c 1-15)
		fi
##Define y (c) coefficient!
		cz1=$( echo "scale=30 ; $cz + $len_z " |bc |cut -c 1-15)
		cz=$( echo "scale=30 ; $cz1"  | bc | cut -c 1-15 )
		czz=$( echo "scale=30 ; $vec_cz - $cz" |bc |cut -c 1-15)
##Print the result!
		if [ $a -le  $num_fix_lay ] ;then
			echo $ca"         "$cb"         "$cz" F F F " >> ./POSCAR
			# echo $caa"         "$cbb"         "$czz" F F F" >> ./POSCAR
			else
			echo $ca"         "$cb"         "$cz" T T T" >> ./POSCAR
			# echo $caa"         "$cbb"         "$czz" T T T" >> ./POSCAR
		fi
	done

}



function poscar_fcc_100_slab_primi() {
	##usage: poscar_fcc_100_slab_primi 'vaccume distance' 'number of atom'
	##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3 ; seq_num=2
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 100 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con / sqrt(2) " | bc |cut -c 1-15)
	vec_by=$( echo "$vec_ax")
	vec_cz=$( echo "scale=30 ; $lat_con / 2 * ( $num_atom - 1 ) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_ax) / 2" |bc )
	len_y=$( echo "scale=30 ; $vec_by / 2" |bc )
	len_z=$( echo "scale=30 ; $lat_con / 2 " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	echo >> ./POSCAR
}

function poscar_fcc_110_slab_primi() {
	##usage: poscar_fcc_110_slab_primi 'vaccume distance'  'number of atom'
		##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  seq_num=2
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim

	echo $material_name 110 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
	vec_ax=$( echo "$lat_con" |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $lat_con / sqrt(2) " | bc |cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con / sqrt( 2 ) / 2 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_ax) / 2" |bc )
	len_y=$( echo "scale=30 ; $vec_by / 2" |bc )
	len_z=$( echo "scale=30 ; $lat_con  / sqrt( 2 ) / 2  " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	echo >> ./POSCAR

}

function poscar_fcc_111_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  seq_num=3;  
	if [[ $4 ]]; then
		strain_A=$4;  
		if [[ $5 ]]; then
		 	tau=$5
		fi 
	else
		strain_A=0;   tau=0
	fi
	strain_l=$( echo "scale=30 ; sqrt( $strain_A + 1 ) - 1 " | bc |cut -c 1-30)
	# lz=$4
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
    echo $material_name 111 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR
	vec_ax=$( echo "scale=30 ; $lat_con / sqrt(2) * ( 1 + $strain_l ) " | bc |cut -c 1-15)
	vec_bx=$( echo "scale=30 ; $vec_ax * 0.5 * ( 1 + $strain_l ) " | bc |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $vec_ax * 0.5 * sqrt(3) * ( 1 + $strain_l ) " | bc | cut -c 1-15)
	# vec_cz=$( echo "scale=30 ; $lz " |bc |cut -c 1-15)
	# vec_cz1=$( echo "scale=30 ; $lat_con * sqrt(3) / 3 * ( 1 - $strain_A * $tau ) " |bc |cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con * sqrt(3) / 3 * ( 1 - $strain_A * $tau ) * ($num_atom - 1) + $1 " |bc |cut -c 1-15)
	# echo $strain_l $vec_ax $vec_bx $vec_by $vec_cz1

	poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_bx + $vec_ax) / 3" |bc )
	len_y=$( echo "scale=30 ; $vec_by / 3" |bc )
	len_z=$( echo "scale=30 ; $lat_con * sqrt(3) / 3 * ( 1 - $strain_A * $tau ) " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num

	echo >> ./POSCAR
}

function poscar_fcc_210_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'
	vaccume=$1; 	num_atom=$2; 	num_rlx=$3 ; seq_num=10
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim

	echo $material_name 210 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(6) * 0.5 " | bc |cut -c 1-15) #done
	vec_bx=$( echo "scale=30 ; $vec_ax * 2 / 3 " | bc |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $vec_ax * sqrt(5) / 3 " | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con * sqrt(5) / 10 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_bx + $vec_ax ) * 0.7" |bc )
	len_y=$( echo "scale=30 ; ($vec_by ) * 0.7" |bc )
	len_z=$( echo "scale=30 ; $lat_con * sqrt(5) / 10 " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num

	echo >> ./POSCAR
}



function poscar_fcc_211_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  seq_num=6
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 211 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(2) * 0.5 " | bc |cut -c 1-15) #done
	vec_by=$( echo "scale=30 ; $lat_con * sqrt(3) " | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con * sqrt(6) * 0.5 / 6 * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax 0.0000000000 $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_ax ) * 0.5 " |bc )
	len_y=$( echo "scale=30 ; ($vec_by ) * 2 / 3" |bc )
	len_z=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 / 6  " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	echo >> ./POSCAR

}


function poscar_fcc_311_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  seq_num=11
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 311 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 " | bc |cut -c 1-15) #done
	vec_bx=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 * 5 / 6 " | bc |cut -c 1-15)  ## Cos (theta) =5/6
	vec_by=$( echo "scale=30 ; $lat_con * sqrt(6) / 2 * sqrt(11) / 6  " | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con / sqrt(11)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_ax + $vec_bx) *  0.2727 " |bc )
	len_y=$( echo "scale=30 ; $vec_by * 0.2727 " |bc )
	len_z=$( echo "scale=30 ; $lat_con / sqrt(11) " |bc )
	cb=0; 	ca=0;	cz=0

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num

	echo >> ./POSCAR
}


function poscar_fcc_331_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  seq_num=19 #seq_num=11
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 331 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 " | bc |cut -c 1-15) #done
	vec_bx=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 * 0.9 " | bc |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $lat_con * sqrt(10) / 2 * sqrt(0.19) " | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con / sqrt(19)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; ($vec_ax + $vec_bx) * ( 1 - sqrt(10) / 10 ) " |bc )
	len_y=$( echo "scale=30 ; $vec_by * ( 1 - sqrt(10) / 10 )  " |bc )
	len_z=$( echo "scale=30 ; $lat_con / sqrt(19) " |bc )

 	scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	echo >> ./POSCAR
}
### Editting unit cell
function poscar_fcc_557_slab_primi(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  asym=$4; seq_num=100
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 557 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(2) / 2 " | bc |cut -c 1-15) #done
	vec_bx=$( echo "scale=30 ; $lat_con * sqrt(2) * 5 / 4 " | bc |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $lat_con * sqrt(22) * 3  / 4" | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con / (3.31662479036*3)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	# poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; $lat_con * 3 / sqrt(2) " |bc )
	len_y=$( echo "scale=30 ; $lat_con * sqrt(127314) / 99 * sqrt(1-(297)^2/(127314*2))  " |bc )
	# len_x=$( echo "scale=30 ; ($vec_ax)* (99-17) / 99  + ($vec_bx) * (99-17) / 99  " |bc )
	# len_y=$( echo "scale=30 ; $vec_by * (99-17) / 99   " |bc )
	len_z=$( echo "scale=30 ; $lat_con / sqrt(5^2+5^2+7^2)  " |bc )

	if [[ $asym == 1 ]]; then
		#this is asymmeteric slab making
		# num_atom2=$((2 * $2 ))
		scale=$num_atom
		vec_cz=$( echo "scale=30 ; $lat_con / (3.31662479036*3)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
		poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
		scaling_making_asym_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	else
		#this is symmeteric slab making
		poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
		scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	fi

	echo >> ./POSCAR
}
### Editting unit cell
function poscar_fcc_557_slab_primi_re(){
        ##usage: poscar_fcc_111_slab_primi 'vaccume distance' 'number of atom'
        ##please use this function with primitive FCC cell caltulation result 'CONTCAR'

	vaccume=$1; 	num_atom=$2; 	num_rlx=$3;  asym=$4; seq_num=100
	poscar_head $num_atom $num_rlx; poscar_lattic_from_contcar_prim
	echo $material_name 557 $num_atom AL > ./POSCAR; echo 1.0 >> ./POSCAR

	vec_ax=$( echo "scale=30 ; $lat_con * sqrt(2) / 2 " | bc |cut -c 1-15) #done
	vec_bx=$( echo "scale=30 ; $lat_con * sqrt(2) * 5 / 4 " | bc |cut -c 1-15)
	vec_by=$( echo "scale=30 ; $lat_con * sqrt(22) * 3  / 4" | bc | cut -c 1-15)
	vec_cz=$( echo "scale=30 ; $lat_con / (3.31662479036*3)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
	# poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz

	test_x=$( echo "scale=30 ; $vec_bx + $vec_ax " |bc )
	test_y=$( echo "scale=30 ; $vec_by " |bc )
	len_x=$( echo "scale=30 ; $lat_con * 3 / sqrt(2) " |bc )
	len_y=$( echo "scale=30 ; $lat_con * sqrt(127314) / 99 * sqrt(1-(297)^2/(127314*2))  " |bc )
	# len_x=$( echo "scale=30 ; ($vec_ax)* (99-17) / 99  + ($vec_bx) * (99-17) / 99  " |bc )
	# len_y=$( echo "scale=30 ; $vec_by * (99-17) / 99   " |bc )
	len_z=$( echo "scale=30 ; $lat_con / sqrt(5^2+5^2+7^2)  " |bc )

	if [[ $asym == 1 ]]; then
		#this is asymmeteric slab making
		# num_atom2=$((2 * $2 ))
		scale=$num_atom
		vec_cz=$( echo "scale=30 ; $lat_con / (3.31662479036*3)  * ($num_atom - 1) + $vaccume " |bc |cut -c 1-15)
		poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
		scaling_making_asym_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	else
		#this is symmeteric slab making
		poscar_slab_vector_selec $vec_ax $vec_bx $vec_by $vec_cz
		scaling_making_slab $test_x $test_y $len_x $len_y $len_z $num_fix_lay $seq_num
	fi

	echo >> ./POSCAR
}
# vac=16
# for a in 11 21 37
# do
# 	num=$a
# 	rlx=$num
# 	for b in -0.04 -0.02 0 0.02 0.04
# 	do
# 		str=$b
# 		tau_hkl=0.45
# 		poscar_fcc_111_slab_primi $vac $num $rlx $str $tau_hkl ; mv POSCAR FCC111_"$str"_"$num"AL.vasp
# 	done
# done


# z_eq=2.29
# poscar_fcc_111_slab_primi $vac $num	$rlx $z_eq 	; mv POSCAR FCC111_$num.vasp
# poscar_fcc_111_slab_primi $vac $num	$rlx	; mv POSCAR FCC111_$num.vasp


# poscar_fcc_557_slab_primi $vac $num	$rlx 1 
poscar_fcc_111_slab_primi $vac $num	$rlx	; mv POSCAR FCC111_$num.vasp
# poscar_fcc_110_slab_primi $vac $num	$rlx	; mv POSCAR FCC110_$num.vasp
# poscar_fcc_100_slab_primi $vac $num	$rlx	; mv POSCAR FCC100_$num.vasp
# poscar_fcc_210_slab_primi $vac $num	$rlx	; mv POSCAR FCC210_$num.vasp
# poscar_fcc_211_slab_primi $vac $num	$rlx	; mv POSCAR FCC211_$num.vasp
# poscar_fcc_311_slab_primi $vac $num	$rlx	; mv POSCAR FCC311_$num.vasp
# poscar_fcc_331_slab_primi $vac $num	$rlx	; mv POSCAR FCC331_$num.vasp
