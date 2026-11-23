#!/bin/bash
#
##	submit with sBatch, monitor with squeue -u jh4, cancel with scancel id#
##	for array job sBatch -a 1-10 OR -a 3,9
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=DarkLit
#SBATCH --output='./DarkLit03/out/DarkLit_%a.out'
module purge
module load matlab/2018a
cat MVARiDarkLit03.m | srun matlab -nodisplay

