#!/bin/bash
#
##	submit with sBatch, monitor with squeue -u jh4, cancel with scancel id#
##	for array job sBatch -a 1-10
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=DarkLit
##SBATCH --mail-type=END
##SBATCH --mail-user=bob.smith@nyu.edu
#SBATCH --output='./DarkLit02/out/DarkLit_%a.out'

module purge
module load matlab/2018a

##cd /scratch/$USER/mytest1
cat MVARiDarkLit02.m | srun matlab -nodisplay

