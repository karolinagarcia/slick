#!/bin/bash

srun slick_init_jobscript.sh $1
srun slick_run_jobscript.sh
echo "Slick has finished"