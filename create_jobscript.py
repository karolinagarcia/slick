def create_jobscript(param_filename, max_lines, SBATCH_args={}):
    default_SBATCH_args = {
        "qos": "narayanan-b",
        "nodes": "1",
        "tasks-per-node": "1",
        "cpus-per-task": "1",
        "mem-per-cpu": "8gb",
        "time": "96:00:00",
    }
    # NOTE: this makes no guarantees about the order in which SBATCH args are laid out
    full_SBATCH_args = {**default_SBATCH_args, **SBATCH_args}
    full_SBATCH_args["array"] = f"1-{max_lines}"

    with open("slick_run_jobscript.sh", "w") as f:
        f.write("#!/bin/bash")
        for arg_name, arg_val in full_SBATCH_args.items():
            f.write(f"#SBATCH --{arg_name}={arg_val}")
        f.write(f"""PRAM=$(sed -n "$SLURM_ARRAY_TASK_ID"p {param_filename})""")
        f.write("""echo $PARAM""")
        f.write("""python slick_run.py --cloudinfo $PARAM""")
