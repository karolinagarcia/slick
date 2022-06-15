#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

python "$parent_path/slick_init.py" $1
if [ $? = 0 ]
then
  srun bash slick_run_jobscript.sh
  echo "Slick has finished"
  exit 0
else
  echo "Slick init failed"
  exit 1
fi
