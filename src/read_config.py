from configparser import ConfigParser, NoSectionError

def parse_parameters(config_file):
    config = ConfigParser()
    config.read(config_file)
    config_result = {}
    config_result["boxsize"] = config['snap']['boxsize']
    config_result["ytfilename"] = config['snap']['ytfilename']
    config_result["caesarfilename"] = config['snap']['caesarfilename']
    config_result["mode"] = config['sample']['mode']
    if config_result["mode"] == "randomize":
        config_result["n_galaxies_sample"] = config['sample']['n_galaxies_sample']
        config_result["min_mass"] = config['sample']['min_mass']
        config_result["max_mass"] = config['sample']['max_mass']
    try:
        config_result["sbatch"] = {k: v for k, v in config.items("sbatch")}
    except NoSectionError:
        config_result["sbatch"] = {}
    config_result["output_dir"] = config.get("run", "output_dir", fallback="Output_Files")
    config_result["basictable_dir"] = config.get("run", "basictable_dir", fallback="Output_Files")
    config_result["skip_lumcalc"] = config.getboolean("run", "skip_lumcalc", fallback=False)
    config_result["skip_basictable"] = config.getboolean("run", "skip_basictable", fallback=False)
    config_result["overwrite"] = config.getboolean("run", "overwrite", fallback=False)
    return config_result
