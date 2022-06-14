#!/bin/bash
#JOBID=$(sbatch slick1_jobscript.sh | awk '{print $4}')
#JOBID=$(sbatch -d afterok:$JOBID slick2_jobscript.sh | awk '{print $4}')

sbatch --dependency=singleton --job-name=slick slick1_jobscript.sh
sbatch --dependency=singleton --job-name=slick slick2_jobscript.sh