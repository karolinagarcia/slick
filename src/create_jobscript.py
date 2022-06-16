import pathlib

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
        f.write("#!/bin/bash\n")
        f.writelines(f"#SBATCH --{arg_name}={arg_val}\n" for arg_name, arg_val in full_SBATCH_args.items())
        f.write(f"python {pathlib.Path(__file__).parent.resolve()}/slick_run.py --cloudinfofile {param_filename} --cloudinfoline $SLURM_ARRAY_TASK_ID\n")
