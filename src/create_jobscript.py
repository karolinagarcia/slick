import pathlib

def create_jobscript(param_filename, max_lines, param_file, SBATCH_args={}):
    default_SBATCH_args = {
        "nodes": "1",
        "tasks-per-node": "1",
        "cpus-per-task": "1",
        "mem-per-cpu": "8gb",
        "time": "48:00:00",
        "qos": "narayanan-b",
        "output": "/dev/null"
    }
    # NOTE: this makes no guarantees about the order in which SBATCH args are laid out
    full_SBATCH_args = {**default_SBATCH_args, **SBATCH_args}
    full_SBATCH_args["array"] = f"1-{max_lines}%92"

    with open("slick_run_jobscript.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.writelines(f"#SBATCH --{arg_name}={arg_val}\n" for arg_name, arg_val in full_SBATCH_args.items())
        f.write(f"python {pathlib.Path(__file__).parent.resolve()}/slick_run.py --parameters {param_file} --cloudinfofile {param_filename} --cloudinfoline $SLURM_ARRAY_TASK_ID\n")
