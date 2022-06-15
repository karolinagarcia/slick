def create_jobscript(out_dir, *, SBATCH_args={}):
    default_SBATCH_args = {
        "qos": "narayanan-b",
        "nodes": "1",
        "tasks-per-node": "1",
        "cpus-per-task": "1",
        "mem-per-cpu": "8gb",
        "time": "96:00:00",
        # NOTE: array param depends on # of clouds
        # "array": f"1-{max_lines}"
    }
    # NOTE: this makes no guarantees about the order in which SBATCH args are laid out
    full_SBATCH_args = {**default_SBATCH_args, **SBATCH_args}

    with open(f"{out_dir}/slick_run_jobscript.sh", "w") as f:
        f.write("#!/bin/bash")
        for arg_name, arg_val in full_SBATCH_args.items():
            f.write(f"#SBATCH --{arg_name}={arg_val}")
        f.write("""PRAM=$(sed -n "$SLURM_ARRAY_TASK_ID"p Clouds_per_Core_m400_z=0.0_total.txt)""")
        f.write("""echo $PARAM""")
        f.write("""python slick_run.py --cloudinfo $PARAM""")
