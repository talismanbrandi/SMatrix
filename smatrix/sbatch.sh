#!/bin/bash
#SBATCH --nodes=1                       # number of nodes requested
#SBATCH --ntasks=20                     # number of tasks (default: 1)
#SBATCH --partition=all                 # partition to run in (all or maxwell)
#SBATCH --job-name=NAME                 # job name
#SBATCH --output=logs/NAME-%N-%j.out    # output file name
#SBATCH --error=logs/NAME-%N-%j.err     # error file name
#SBATCH --time=1:00:00                  # runtime requested
#SBATCH --mail-user=ayan.paul@desy.de   # notification email
#SBATCH --mail-type=FAIL                # notification type
#SBATCH --export=ALL
#SBATCH --array 1-2048
export LD_PRELOAD=""

# run the application:
number=$((${SLURM_ARRAY_TASK_ID}-1))
cd configuration_${number}
mpiexec ./analysis config/myModel.conf config/MonteCarlo.conf
cd ..
