#!/bin/sh

#SBATCH --job-name=JOBNAME
#SBATCH --output=LOGFILENAME
#SBATCH --error=ERRORFILENAME
#SBATCH --partition=htc
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3000
#SBATCH --time=03:00:00
#SBATCH --licenses=sps
#SBATCH --account=hyperk

start=$(date +%s)

TASKS

end=$(date +%s)
echo "Time JOBNAME: $((($end-$start)/60)) mins"
