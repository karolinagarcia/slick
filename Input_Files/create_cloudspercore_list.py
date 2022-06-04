import yt
import caesar
import random
import numpy as np
from tqdm import tqdm

def create_param_file(boxsize,ytfilename,caesarfilename,mode,n_galaxies_sample,min_mass,max_mass):
    yt_snap305 = yt.load(ytfilename)
    yt_snap305_data = yt_snap305.all_data()
    obj = caesar.load(caesarfilename)

    if mode == 'randomize':
        galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies if (gal.masses['total']>min_mass and gal.masses['total']<=max_mass)]
        galaxies_list = random.sample(galaxies_list, n_galaxies_sample)
    else: galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies]

    n_of_gals = len(galaxies_list)
    n_of_clouds = sum([len(galaxies_list[i][1]) for i in np.arange(n_of_gals)])
    n_clouds_per_line = n_of_clouds/700

    #with open('parameters_'+str(n_galaxies_sample)+'g_'+str(n_clouds_per_line)+'cperline_mass'+str(min_mass)+'_to_'+str(max_mass)+'.txt', 'w') as f:
    paramfilename = '../Input_Files/parameters_boxlen='+str(boxsize)+'Mpc-h_z='+str(round(yt_snap305.parameters['Redshift'],3))+'_'+mode+'.txt'
    with open(paramfilename, 'w') as f:
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
        f.close()

    return paramfilename
