#!/bin/sh

# Set the global variables:
while getopts F:g:R:r:t:n:N:T:S: flag
do
	case "${flag}" in
		F) folder=${OPTARG};;  # F - folder name
		g) geo=${OPTARG};;     # g - geometry ["cave", "chimney"]
		R) R=${OPTARG};;       # R - the initial radius of the geometry.
		r) r=${OPTARG};;       # r - the plume radius / radius of affect w-o radial dependency
		t) dt=${OPTARG};;      # t - the timestep of each prognostic iteration [hours]
		n) nt=${OPTARG};;      # n - the number of prognostic iterations
		N) N=${OPTARG};;       # N - number of iterations
		T) T=${OPTARG};;       # T - The total run time allowed [s]
		S) T_i=${OPTARG};;	   # S - Start time
	esac
done

# Extract the numbers from the folder names.
IFS='_' read -ra array <<< "$folder"
IFS='H' read -ra thickness <<< "${array[0]}"
H=${thickness[1]}
IFS='I' read -ra inclination <<< "${array[1]}"
I=${inclination[1]}
IFS='Q' read -ra heat_flux <<< "${array[2]}"
Q=${heat_flux[1]}
IFS='M' read -ra mode <<< "${array[3]}"
M=${mode[1]}

# Convert the length of each prognostic time step from hours to years.
dt=`python -c "print(${dt}/24/365.25)"`  # We need python to handle the floating point number.

# Change the total geothermal heat flux from kW to W.
Q=$((1000*$Q))

# The variables needed to infill the DEM parameters in the .sif files.
dx=`python -c 'print(0.5)'`
x0=`python -c 'print(-1.0)'`
y0=`python -c 'print(-1.0)'`
lx=`python -c "print('%.1f' % (2*(${H}+1)))"`
ly=`python -c "print('%.1f' % (2*(${H}+1)))"`
Nx=`python -c "print('%d' % (2*(${H}+1)/${dx} + 1))"`
Ny=`python -c "print('%d' % (2*(${H}+1)/${dx} + 1))"`

# Parent directory (pdir) and current working directory (cwd).
pdir=$(pwd)
cwd="${pdir}/${folder}"
cd $cwd

# Directories
DIR_geo_results="${cwd}/${geo}_results"
DIR_geo="${cwd}/${geo}"
DIR_temp="${cwd}/temp"

# File paths
f_mesh="${cwd}/${geo}.msh"
f_csv="${cwd}/${geo}.csv"
f_diagnostic_sif_fillin="${pdir}/main/diagnostic.sif"
f_prognostic_sif_fillin="${pdir}/main/prognostic.sif"
f_diagnostic_sif="${cwd}/diagnostic.sif"
f_prognostic_sif="${cwd}/prognostic.sif"
f_diagnostic_result="${DIR_geo}/${geo}_diagnostic_t0001.vtu"
f_prognostic_result="${DIR_geo}/${geo}_prognostic_t000${nt}.vtu"
f_mesh_void="${DIR_temp}/surf101.stl"

# The python files.
p_vtu2msh="${pdir}/main/vtu2msh.py"
p_remesh="${pdir}/main/remesh.py"



##############################
###### --- Functions --- #####
##############################

function fill_in_sif {
	# Fill in the appropreate values in the solver input files.
	# $1 - input file path
	# $2 - output file path
	sed -e "s/#geo#/${geo}/g"\
		-e "s/#nt#/${nt}/g"\
		-e "s/#time#/${dt}/g"\
		-e "s/#I#/${I}/g"\
		-e "s/#x0#/${x0}/g"\
		-e "s/#y0#/${y0}/g"\
		-e "s/#lx#/${lx}/g"\
		-e "s/#ly#/${ly}/g"\
		-e "s/#Nx#/${Nx}/g"\
		-e "s/#Ny#/${Ny}/g"	$1 > $2
}


function mesh_resolution {
	# Set the mesh resolutions.
	if [ $H -le 50 ]
	then
		minRes=3
		maxRes=10
	elif [ $H -le 100 ]
	then
		minRes=3
		maxRes=15
	elif [ $H -le 150 ]
	then
		minRes=3
		maxRes=20
	elif [ $H -le 200 ]
	then
		minRes=3
		maxRes=30
	else
		minRes=3
		maxRes=40
	fi
}


function diagnostic_run {
	local start_time=`date +%s`
	echo "- Diagnostic run"
	
	# Elmer/Ice diagnostic run.
	ElmerGrid 14 2 "${geo}.msh" -autoclean >/dev/null
	ElmerSolver diagnostic.sif >/dev/null
	
	local end_time=`date +%s`
	echo " -- Diagnostic runtime:" $((end_time-start_time))
}


function prognostic_run {
	local start_time=`date +%s`
	echo "- Prognostic run"
	
	local prognostic_condition=1
	
	while [ $prognostic_condition -eq 1 ]
	do
		# Elmer/Ice prognostic run.
		ElmerSolver prognostic.sif >/dev/null
		
		# Check if the prognostic run was successful.
		prognostic_check $f_prognostic_result
		PrognosticCheck=$?
		
		if test ! -e $f_prognostic_result
		then
			echo "The prognostic run failed for some reason."
			exit 1
		elif [ $PrognosticCheck -eq 0 ]
		then
			# Shorten the time step.
			echo " -- Timestep possibly too large. Trying 3/4*dt..."
			dt=`python -c "print(3*${dt}/4)"`
			
			if [ `python -c "print(1 if ${dt} < 1/2/24/365.25 else 0)"` -eq 1 ]
			then
				" -- Timestep is not the problem."
				prognostic_condition=0
				condition=0
			fi
			
			fill_in_sif	$f_prognostic_sif_fillin $f_prognostic_sif
		else
			local prognostic_condition=0
		fi
	done
	
	local end_time=`date +%s`
	echo " -- Prognostic runtime:" $((end_time-start_time))
}


function vtu2msh {
	local start_time=`date +%s`
	
	# Remesh the input diagnostic result ($1).
	python $p_vtu2msh -g $geo -H $H -I $I -f $1 -m $minRes -M $maxRes
	
	local end_time=`date +%s`
	echo " -- Remesh runtime:" $((end_time-start_time))
}


function remesh {
	local start_time=`date +%s`
	
	# Deploy the remeshing scheme on the desired inpute prognostic results ($1).
	python $p_remesh -g $geo -H $H -I $I -r $r -f $1 -m $minRes -M $maxRes -Q $Q -QM $M -dt $dt -nt $nt
	
	local end_time=`date +%s`
	echo " -- Remesh runtime:" $((end_time-start_time))
	
	# Check if a mesh was created during the remeshing.
	if test ! -f $f_mesh
	then
		condition=1
		echo "Mesh was not created."
	fi
}


function prognostic_check {
	# Check if numerical instabilities have resulted in weirdly shaped mesh.
	local min_val=`python -c "import pyvista as pv; vtu=pv.read('$1'); print(vtu.points[:, 2].min())"`
	local max_val=`python -c "import pyvista as pv; vtu=pv.read('$1'); print(vtu.points[:, 2].max())"`
	local max_diff=`python -c "import math; print((1.05+2*math.tan(${I}*math.pi/180))*${H})"`
	local diff_condition=`python -c "print(1 if (${max_val} - ${min_val}) < ${max_diff} else 0)"`
	
	return $diff_condition
}


function condition_check {
	# Change the condition if any parameters break their limits.
	local condition_small=`python -c "import numpy;\
						   print(1 if numpy.array([${h_z}, ${h_a}, ${h_c}]).min() < 2 else 0)"`
	local condition_large=`python -c "print(1 if ${d_i} < ${H}/4 else 0)"`
	local condition_width=`python -c "print(1 if ${x_v}/2 > 3*${H}/4 or ${y_v}/2 > 3*${H}/4 else 0)"`
	local T_now=`date +%s`
	
	if [ $i -ge $N ]
	then
		condition=1
		echo "Maximum number of iterations reached (${N})."
	elif [ $condition_small -eq 1 ]
	then
		echo ${h_z} ${h_a} ${h_c}
		condition=1
		echo "Void too small / closed."
	#elif [ $condition_large -eq 1 ]
	#then
	#	condition=1
	#	echo "Void has reached within a quarter of the ice surface."
	elif [ $condition_width -eq 1 ]
	then
		condition=1
		echo "Void footprint has grown too wide."
	elif [ $(($T_now-$T_i)) -ge $T ]
	then
		condition=1
		echo "Run-time limit reached."
	fi
}


function print_head {
	local message="Step $1:"
	local border=""
	for j in $(seq 1 ${#message}); do border+="="; done

	echo $border
	echo $message
	echo $border
}


function num2str {
	# Create a number string that handles four significant digits: 
	# 	i.e. 1 -> 0001, 12 -> 0012, 123 -> 0123 and 1234 -> 1234.
	if [ $1 -lt 10 ]
	then
		num="000$1"
	elif [ $1 -gt 9 ] && [ $1 -lt 100 ]
	then
		num="00$1"
	elif [ $1 -gt 99 ] && [ $1 -lt 1000 ]
	then
		num="0$1"
	else
		num="$1"
	fi
}



#####################################################################################################
# Step n: Extract the mesh from the last diagnostic run and run both diagnostic and prognostic runs #
#####################################################################################################

# Find the last diagnostic vtu file.
VTUs=`find "${DIR_geo_results}" -type f -name '*.vtu*'`
IFS=$'\n' vtus_sorted=($(sort <<<"${VTUs[*]}")); unset IFS
f_diagnostic_n=${vtus_sorted[-1]}

# Extract the iteration number from from the .vtu file name.
IFS='_' read -ra array <<< $f_diagnostic_n; unset IFS
IFS='.' read -ra suffix <<< "${array[-1]}"; unset IFS
IFS='t' read -ra nn <<< "${suffix[0]}"; unset IFS
IFS='0' read -ra nn <<< "${nn[-1]}"; unset IFS
nn=${nn[-1]}
NN=0
for (( i=0; i<${#nn}; i++ ))
do
	NN=$(($NN + 10**(${#nn}-1-i) * ${nn:i:1}))
done

num2str $NN

# Determine the appropreate mesh resolutions
mesh_resolution

# Print prelogue.
print_head $nn
echo `date +%s`
echo ""
echo "Thickness: ${H} m"
echo "Bed slope: ${I} degrees"
Q_str=`python -c "print('%.2f' % (${Q}*1e-6))"`
echo "Heat flux: ${Q_str} MW"
echo "mode: ${M}"
echo " "


# Copy and modify the .sif files to contain the correct geometry name.
fill_in_sif $f_diagnostic_sif_fillin $f_diagnostic_sif
fill_in_sif $f_prognostic_sif_fillin $f_prognostic_sif


# Get the mesh from the last diagnostic run.
vtu2msh $f_diagnostic_n

# Diagnostic run
diagnostic_run

# Copy the diagnostic solution to a folder that will keep all the results.
cp $f_diagnostic_result $f_diagnostic_n

# Prognostic run
prognostic_run



############################################################################
# Step n+: Create a mesh from the previous output and run d-c and p-c runs #
############################################################################

# Now iterate over the desired timespan.
i=$(($nn+1))
condition=0
while  [ $condition -eq 0 ]
do
	# Print prelogue.
	print_head $i
	echo `date +%s`
	
	# Deploy the remeshing scheme:
	remesh $f_prognostic_result
	
	# Delete the folders (and their files) created during the remeshing and E/I operations.
	#rm -r $DIR_geo
	
	# Diagnostic run
	if [ $condition -eq 0 ]
	then
		diagnostic_run
	
		# Copy the results to the apppropreate folder
		num2str $i
		f_diagnostic_i="${DIR_geo_results}/${geo}_diagnostic_t${num}.vtu"
		cp $f_diagnostic_result $f_diagnostic_i
	
		# Read the new variables from the .csv file
		IFS=, read -r V_v A_v h_z h_a h_c x_v y_v q V_i d_i H_i t < <(tail -n1 $f_csv)
	fi
	
	# Check if any parameters have reached the limits.
	condition_check
	
	# Prognostic run (do not run in the final time step).
	if [ $condition -eq 0 ]
	then
		#Prognostic run.
		prognostic_run
	fi
	
	# Remove the mesh if there is to be another iteration.
	if [ $condition -eq 0 ]
	then
		# Delete the mesh.
		rm -r $f_mesh
	fi
	
	((i++))
done

cd ..
	
