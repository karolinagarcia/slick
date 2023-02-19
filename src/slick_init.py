from sys import argv
from os import mkdir

from create_basic_characteristics_table import create_basic_table
from create_cloudspercore_list import create_cloudspercore_list
from create_jobscript import create_jobscript
from read_config import parse_parameters

def main():
    if len(argv) > 1:
        config_file = argv[1]
    else:
        print("no parameter file given.")
        # TODO: print some sort of usage message here?
        exit(1)

    config = parse_parameters(config_file)
    
    #mkdir(config["output_dir"])

    # Creates table with basic characteristics for all the clouds
    if not config["skip_basictable"]:
        create_basic_table(config)

    if not config["skip_lumcalc"]:
        # Creates table with basic characteristics for all the clouds (if mode='total'), or for a sample of them (mode = 'randomize')
        param_filename, max_lines = create_cloudspercore_list(config)

        # Creates slick_run_jobscript.sh
        create_jobscript(param_filename, max_lines, config_file, config["sbatch"])

if __name__ == "__main__":
    main()