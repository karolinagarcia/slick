#!/bin/bash
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

python "$parent_path/../src/slick_init.py" $1
if [ $? = 0 ]
then
  sbatch slick_run_jobscript.sh $1
  echo "Slick array queued"
  exit 0
else
  echo "Slick init failed"
  exit 1
fi
