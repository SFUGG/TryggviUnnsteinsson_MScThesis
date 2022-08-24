#!/bin/sh

# The maximum number of iterations (N=600 gives roughly six months).
N=#N#
# The void geometry type.
geo="#geo#"
# The initial geometry radius.
R=#R#
# The plume radius.
r=#r#
# The length of each prognostic time step [hours].
dt=#dt#
# Number of prognostic timesteps.
nt=#nt#

# The maximum run time (in seconds).
hrs=#hrs#
runtime=$(($hrs*60*60))

# Get a list of all the sub-folders and delete the last entry (main folder).
folder_list=$(while IFS= read -r DIR; do echo "${DIR}"; done < <(find . -maxdepth 1 -mindepth 1 -type d -printf "%P\n"))
folder_list=(${folder_list[@]/main})

# Iterate through all the folders.
for DIR in ${folder_list[*]}
do
	start=`date +%s`
	bash main/glaciovolcanic_void.sh -F $DIR -g $geo -R $R -r $r -t $dt -n $nt -N $N -T $runtime -S `date +%s`
	end=`date +%s`
	sec_TOT=$(($end - $start))
	sec=$((sec_TOT % 60))
	min=$(((sec_TOT - sec)/60 % 60))
	hrs=$((((sec_TOT - sec)/60 - min)/60))
	echo "Runtime ${folder}: ${hrs}:${min}:${sec}"
done
