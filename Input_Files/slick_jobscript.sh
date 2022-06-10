#!/bin/bash
#
#SBATCH --qos=narayanan-b
#SBATCH --output=../Output_Tables/SLURM_20220603
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --array=1-2942

PRAM=$(sed -n "$SLURM_ARRAY_TASK_ID"p Clouds_per_Core_m100_z=0.993_total.txt)
echo $PRAM

python slick.py --cloudinfo $PRAM
