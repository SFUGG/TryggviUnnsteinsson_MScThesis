#!/bin/sh
#SBATCH --time #hrs#:00:00
#SBATCH --mem=4000M
#SBATCH --account=#account#
#SBATCH --array=1-#total#

# Load Elmer
module load gcc/9.3.0
module load elmerfem/9.0

# Load python
module load python/3.6

# Load python virtual environment
source $HOME/elmericepython/bin/activate

echo "Starting task $SLURM_ARRAY_TASK_ID"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" folder_list.txt)

# For timing.
start=`date +%s`

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
hrs=#hrs_1#
runtime=$(($hrs*60*60))

bash main/glaciovolcanic_void_restart.sh -F $DIR -g $geo -R $R -r $r -t $dt -n $nt -N $N -T $runtime -S `date +%s`

deactivate

# For timing.
end=`date +%s`
sec_TOT=$(($end - $start))
sec=$((sec_TOT % 60))
min=$(((sec_TOT - sec)/60 % 60))
hrs=$((((sec_TOT - sec)/60 - min)/60))
echo $start
echo $end
echo "Runtime ${DIR}: ${hrs}:${min}:${sec}"
