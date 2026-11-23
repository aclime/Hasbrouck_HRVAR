#!/bin/bash
#
##	submit with sBatch, monitor with squeue -u jh4, cancel with scancel id#
##	for array job sBatch -a 1-10
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=16GB

#SBATCH --job-name=EventTime
#SBATCH --output='./EventTime/out/et_%a.out'
module purge
module load matlab/2018a
cat MVARiPartSIP01ET.m | srun matlab -nodisplay
cat  MVARiLexExl01ET.m | srun matlab -nodisplay
cat MVARiDarkLit01ET.m | srun matlab -nodisplay

