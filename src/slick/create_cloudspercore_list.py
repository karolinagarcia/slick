import numpy as np
import pandas as pd
import caesar
import random
random.seed(10)
from .create_basic_characteristics_table import get
import h5py
import glob


def create_cloudspercore_list(config, n_clouds_per_core):
    if config["sim"] == 'SIMBA':
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

        param_filename = f"{config['output_dir']}/Clouds_per_Core_SIMBA{config['boxsize']}_{config['output_name']}.txt"
        with open(param_filename, "w") as f:
            num_lines = 1
            gal_clouds = np.concatenate(
                np.array([gal[1] for gal in galaxies_list], dtype=object)
            )
            try:
                pd_clouds_with_lum = pd.read_csv(f"{config['output_dir']}/lim_df.csv")
                ids_clouds_with_lum = pd_clouds_with_lum["Cloud_ID"]
                #all_clouds = pd.read_csv(glob.glob(f"{config['basictable_dir']}/Basic*")[0])
                all_clouds = pd.read_csv(glob.glob(f"{config['basictable_dir']}/Basic_Characteristics_SIMBA*{config['output_name']}.csv")[0])
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

    if config["sim"] == 'TNG':
        snap_prog_url = "http://www.tng-project.org/api/"+config["sim_code"]+"/snapshots/"+config["snap_number"]+"/"

        gal_ids = config["subhalo_number"].strip('[]').split(',')
        gal_ids = [int(num) for num in gal_ids]

        part_id = np.array([])

        for gal_id in gal_ids:
            sub_prog_url = snap_prog_url+"subhalos/"+str(gal_id)+"/"

            cutout_request = {'gas':'ParticleIDs'}
            cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)
            part_id_aux = h5py.File(cutout,'r')['PartType0']['ParticleIDs'][:]
            part_id = np.concatenate((part_id, part_id_aux))
            part_id = part_id.astype(int)

        #param_filename = f'{config["output_dir"]}/Clouds_per_Core_{config["sim_code"]}_snap{config["snap_number"]}_subhalo{config["subhalo_number"]}.txt'
        param_filename = f'{config["output_dir"]}/Clouds_per_Core_{config["sim_code"]}_snap{config["snap_number"]}_{config["output_name"]}.txt'

        with open(param_filename, 'w') as f:
            num_lines = 1
            gal_clouds = part_id
            try:
                pd_clouds_with_lum = pd.read_csv(f'{config["output_dir"]}/lim_df.csv')
                ids_clouds_with_lum = pd_clouds_with_lum['Cloud_ID']
                #all_clouds = pd.read_csv(glob.glob(f'{config["basictable_dir"]}/Basic*')[0])
                all_clouds = pd.read_csv(f'{config["basictable_dir"]}/Basic_Characteristics_'+config["sim_code"]+'_snap'+config["snap_number"]+'_'+config["output_name"]+'.csv')
                all_ids = all_clouds['c_Index']
                gal_clouds = np.array(all_ids[~all_ids.isin(ids_clouds_with_lum)])
            except: pass
            for (i,), cloud in np.ndenumerate(gal_clouds):
                if i % int(n_clouds_per_core) == 0:
                    if i != 0:
                        f.write("\n")
                        num_lines += 1
                f.write(f"{cloud} ")


    elif config["sim"] == 'CAMELS':

        f = h5py.File("/blue/narayanan/karolina.garcia/github/slick_CAMELS/FOF_Subfind/"+config["sim_type"]+"/"+config["sim_code"]+"/"+config["sim_code_2"]+"/fof_subhalo_tab_0"+config["snap_number"]+".hdf5", 'r')
        index       = int(config["subhalo_number"])
        len_sh      = f['Subhalo/SubhaloLen'][:]
        IDs_sh      = f['IDs/ID'][:]
        start       = np.sum(len_sh[:index])
        end         = start+len_sh[index]
        indexes_sh  = IDs_sh[start:end]
        indexes1 = np.intersect1d(IDs_sh, indexes_sh, assume_unique=False, return_indices=True)[1]
        f.close()
        f = h5py.File("/blue/narayanan/karolina.garcia/github/slick_CAMELS/Sims/"+config["sim_type"]+"/"+config["sim_code"]+"/"+config["sim_code_2"]+"/snap_0"+config["snap_number"]+".hdf5", 'r')
        part_id     = f['PartType0/ParticleIDs'][:][indexes1]
        f.close()

        param_filename = f'{config["output_dir"]}/Clouds_per_Core_{config["sim_type"]}_{config["sim_code_2"]}_snap{config["snap_number"]}_{config['output_name']}.txt'

        with open(param_filename, 'w') as f:
            num_lines = 1
            gal_clouds = part_id
            for (i,), cloud in np.ndenumerate(gal_clouds):
                if i % int(n_clouds_per_core) == 0:
                    if i != 0:
                        f.write("\n")
                        num_lines += 1
                f.write(f"{cloud} ")

    return param_filename, num_lines