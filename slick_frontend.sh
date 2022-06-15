#!/bin/bash

INIT_JOB_ID=$(sbatch --parsable slick_init_jobscript.sh)

sbatch --dependency=afterok:$INIT_JOB_ID slick_run_jobscript.sh

echo "slurm jobs queued"