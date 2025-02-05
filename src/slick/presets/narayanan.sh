#!/bin/bash
#SBATCH --qos=narayanan-b
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=1:00:00

module purge
module load gcc/12.2.0 gsl/2.7

slick init parameters.ini
