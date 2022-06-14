#!/bin/bash
#
#SBATCH --qos=narayanan-b
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --array=1-2000

PRAM=$(sed -n "$SLURM_ARRAY_TASK_ID"p Clouds_per_Core_m400_z=0.0_total.txt)
echo $PRAM

python slick2.py --cloudinfo $PRAM