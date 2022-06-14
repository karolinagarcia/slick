#!/bin/bash
#
#SBATCH --qos=narayanan-b
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=4:00:00

python slick1.py