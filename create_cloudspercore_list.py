import numpy as np
import caesar
import yt
from tqdm import tqdm
import random
random.seed(10)

def create_cloudspercore_list(config):
    yt_snap305 = yt.load(config["ytfilename"])
    yt_snap305_data = yt_snap305.all_data()
    obj = caesar.load(config["caesarfilename"])

    if config["mode"] == 'randomize':
        galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies if (gal.masses['total']>config["min_mass"] and gal.masses['total']<=config["max_mass"])]
        galaxies_list = random.sample(galaxies_list, config["n_galaxies_sample"])
    else: galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies]

    n_of_gals = len(galaxies_list)
    n_of_clouds = sum([len(galaxies_list[i][1]) for i in np.arange(n_of_gals)])
    n_clouds_per_line = n_of_clouds/700

    param_filename = f'Output_Files/Clouds_per_Core_m{config["boxsize"]}_z={round(yt_snap305.parameters["Redshift"],3)}_{config["mode"]}.txt'
    with open(param_filename, 'w') as f:
        num_lines = 1
        for g in tqdm(np.arange(len(galaxies_list))):
            gal_clouds = galaxies_list[g][1]
            sfr_gal = np.sum(yt_snap305_data['PartType0','StarFormationRate'][gal_clouds].value)
            if sfr_gal>0:
                f.write(str(galaxies_list[g][0])+' ')
                batch = 0
                for c in np.arange(len(gal_clouds)):
                    batch+=1
                    if batch > n_clouds_per_line:
                        batch = 1
                        f.write('\n')
                        f.write(str(galaxies_list[g][0])+ ' ')
                    f.write(str(galaxies_list[g][1][c])+' ')
                f.write('\n')
                num_lines += 1
        f.close()
        
    return param_filename, num_lines
