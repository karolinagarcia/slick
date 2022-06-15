from sys import argv
from os import chdir
from configparser import ConfigParser

from create_basic_characteristics_table import create_basic_table
from create_cloudspercore_list import create_cloudspercore_list
from create_jobscript import create_jobscript

def main():
    if len(argv) > 1:
        config_file = argv[1]
    else:
        print("no parameter file given.")
        # TODO: print some sort of usage message here?
        exit(1)

    config = parse_parameters(config_file)

    chdir(config["out_dir"])

    # Creates table with basic characteristics for all the clouds
    create_basic_table(config)

    # Creates table with basic characteristics for all the clouds (if mode='total'), or for a sample of them (mode = 'randomize')
    param_filename, max_lines = create_cloudspercore_list(config)

    # Creates slick_run_jobscript.sh
    create_jobscript(param_filename, max_lines)

def parse_parameters(config_file):
    config = ConfigParser()
    config.read(config_file)
    config_result = {}
    try:
        config_result["out_dir"] = config["out_dir"]["out_dir"]
    except KeyError:
        config_result["out_dir"] = "."
    config_result["boxsize"] = config['snap']['boxsize']
    config_result["ytfilename"] = config['snap']['ytfilename']
    config_result["caesarfilename"] = config['snap']['caesarfilename']
    config_result["mode"] = config['sample']['mode']
    if config_result["mode"] == "randomized":
        config_result["n_galaxies_sample"] = config['sample']['n_galaxies_sample']
        config_result["min_mass"] = config['sample']['min_mass']
        config_result["max_mass"] = config['sample']['max_mass']

if __name__ == "__main__":
    main()