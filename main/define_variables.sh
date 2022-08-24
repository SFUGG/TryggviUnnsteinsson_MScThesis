#!/bin/sh

# Fill in the arrays with the desired values for the following variables.

# Glacier thickness [m]
thickness=(50 100 150 200)
# The bed slope angle [degrees]
slope_angle=(0 5 10 15)
# The total geothermal heat flux [kW]
heat_flux=(0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000)
# The specific modes of geothermal melt applied.
mode=(0 1 2 3)

# The initial geometry of the void.
geo="cave"
# The initial geometry radius.
R=10
# The plume radius.
r=5
# The timestep of prognostic iterations (in hours).
dt=6
# Number of prognostic timesteps.
nt=4
# The maximum number of iterations, set to fill 400 days of real time.
N=`python -c "print(int(400*24/(${dt}*4)))"`

# The total runtime allowed for each simulation.
hrs=13
hrs_1=$(($hrs-1))

# Max number of simultaniously running simulations (only applicaple for running on clusters).
max_number_simulations=32

# Current working directory (cwd).
cwd=$(pwd)

# A file containing the names of all the subfolders.
folder_list="${cwd}/folder_list.txt"

# Create the necessary folders for each experiment.
for H in ${thickness[*]}
do
	
	for I in ${slope_angle[*]}
	do
		
		for Q in ${heat_flux[*]}
		do
		
			for M in ${mode[*]}
			do
				name="H$H"
				name+="_I$I"
				name+="_Q$Q"
				name+="_M$M"
				echo $name
				folder_name="$cwd/$name"
				if test ! -e $folder_name
				then
					mkdir $folder_name
				fi
				
				echo $name >> $folder_list
			done
			
		done

	done

done

# The total number of simulations.
total=$((${#thickness[@]}*${#slope_angle[@]}*${#heat_flux[@]}*${#mode[@]}))
# The number of simulations that can run simultaniously.
if [ $total -lt $max_number_simulations ]
then
	max_number_simulations=$total
fi


function fill_in_sh {
	# Fill in the appropreate values in the bash files.
	# $1 - input file path
	# $2 - output file path
	sed -e "s/#total#/${total}/g"\
		-e "s/#fraction#/${max_number_simulations}/g"\
		-e "s/#hrs#/${hrs}/g"\
		-e "s/#hrs_1#/${hrs_1}/g"\
		-e "s/#N#/${N}/g"\
		-e "s/#geo#/${geo}/g"\
		-e "s/#R#/${R}/g"\
		-e "s/#r#/${r}/g"\
		-e "s/#dt#/${dt}/g"\
		-e "s/#nt#/${nt}/g"	$1 > $2
}

# Fill in for the variables in the job file.
job_in="${cwd}/main/elmer_job_fillin.sh"
job_out="${cwd}/main/elmer_job.sh"
fill_in_sh $job_in $job_out

# Fill in for the variables in the restart job file.
restart_in="${cwd}/main/elmer_job_restart_fillin.sh"
restart_out="${cwd}/main/elmer_job_restart.sh"
fill_in_sh $restart_in $restart_out

# Fill in for the variables in the run_iterative file.
run_in="${cwd}/main/run_iterative_fillin.sh"
run_out="${cwd}/main/run_iterative.sh"
fill_in_sh $run_in $run_out
