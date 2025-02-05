import numpy as np
import pandas as pd
import caesar
import yt
from tqdm import tqdm
import random
import glob

random.seed(10)
import ast


def create_cloudspercore_list(config, n_clouds_per_core):
    obj = caesar.load(config["caesarfilename"])

    if config["mode"] == "randomize":
        galaxies_list = [
            (gal.GroupID, gal.glist)
            for gal in obj.galaxies
            if (
                gal.masses["total"] > config["min_mass"]
                and gal.masses["total"] <= config["max_mass"]
            )
        ]
        galaxies_list = random.sample(galaxies_list, int(config["n_galaxies_sample"]))
    # elif config["mode"] == 'centergal':
    # galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies][0]
    elif config["mode"] == "total":
        galaxies_list = [(gal.GroupID, gal.glist) for gal in obj.galaxies]
    elif config["mode"] == "single":
        galaxies_list = [(gal.GroupID, gal.glist) for gal in obj.galaxies]
        galaxies_list = [galaxies_list[i] for i in ast.literal_eval(config["gal_ids"])]

    # n_of_gals = len(galaxies_list)
    # n_of_clouds = sum([len(galaxies_list[i][1]) for i in np.arange(n_of_gals)])

    param_filename = f"{config['output_dir']}/Clouds_per_Core.txt"
    with open(param_filename, "w") as f:
        num_lines = 1
        gal_clouds = np.concatenate(
            np.array([gal[1] for gal in galaxies_list], dtype=object)
        )
        try:
            pd_clouds_with_lum = pd.read_csv(f"{config['output_dir']}/lim_df.csv")
            ids_clouds_with_lum = pd_clouds_with_lum["Cloud_ID"]
            all_clouds = pd.read_csv(glob.glob(f"{config['basictable_dir']}/Basic*")[0])
            all_ids = all_clouds["c_Index"]
            gal_clouds = np.array(all_ids[~all_ids.isin(ids_clouds_with_lum)])
        except:
            pass
        if len(gal_clouds.flatten()) / int(n_clouds_per_core) > 3000:
            n_clouds_per_core = int(len(gal_clouds.flatten()) / 3000)
        for (i,), cloud in np.ndenumerate(gal_clouds.flatten()):
            if i % int(n_clouds_per_core) == 0:
                if i != 0:
                    f.write("\n")
                    num_lines += 1
            f.write(f"{cloud} ")
    return param_filename, num_lines
