#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

srun $parent_path/slick_init_jobscript.sh $1
srun slick_run_jobscript.sh
echo "Slick has finished"