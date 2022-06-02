# SET THE VARIABLES BELOW:

date = '20220512'

boxsize = 400

ytfilename = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5'
#ytfilename = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_212.hdf5'

caesarfilename = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5'
#caesarfilename = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0212_z1.007.hdf5'


#boxsize = 1000

#ytfilename = '/orange/narayanan/karolina.garcia/simba/m100n1024/snaps/snap_m100n1024_105.hdf5.1'
#caesarfilename = '/orange/narayanan/d.zimmerman/simba/m100n1024/caesar_cats/caesar_simba_105.hdf5'


mode = 'total'
#---
#mode = 'randomize'
#n_galaxies_sample = 30
#min_mass = 10**11
#max_mass = 3*10**12
#---

#---------------------------------------#


import numpy as np
import caesar
import yt
import pandas as pd
import seaborn as sns
import sys
import yt.units as u
import argparse
from tqdm import tqdm
import random
random.seed(10)

from create_param_file import create_param_file
from create_basic_characteristics_table import create_basic_table
from limfunctions import densityProfile, submm_luminosity, creating_table


# Creates table with basic characteristics for all the clouds
basicfilename = create_basic_table(boxsize,ytfilename,caesarfilename)

###

# Creates table with basic characteristics for all the clouds (if mode='total'), or for a sample of them (mode = 'randomize')
paramfilename = create_param_file(boxsize,ytfilename,caesarfilename,mode,n_galaxies_sample,min_mass,max_mass)

###

# Creates final luminosity table
df_basic = pd.read_csv('../Output_Tables/'+basicfilename)

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-ci","--cloudinfo", type=int, nargs="+")

args = parser.parse_args()

gal_number = args.cloudinfo[0]
cloud_list = args.cloudinfo[1:]

df = pd.DataFrame({'Galaxy_ID':[], 'Cloud_ID':[], 'Mcloud':[], 'Rcloud':[], 'Metallicity':[], 'RadField':[], 'redshift':[], 'H2_lcii':[], 'CO10':[], 'CO21':[], 'CO32':[], 'CO43':[], 'CO54':[], 'CO10_intTB':[], 'CO21_intTB':[], 'CO32_intTB':[], 'CO43_intTB':[], 'CO54_intTB':[], 'CI10':[], 'CI21':[], 'CO65':[], 'CO76':[], 'CO87':[], 'CO98':[], 'CO65_intTB':[], 'CO76_intTB':[], 'CO87_intTB':[], 'CO98_intTB':[], 'OI1':[], 'OI2':[], 'OI3':[], 'fH2':[]})

#df.to_csv('/orange/narayanan/karolina.garcia/SLURM_'+date+'/lim_df_'+date+'.csv', index = False)
df.to_csv('../Output_Tables/SLURM_'+date+'/lim_df_'+date+'.csv', index = False)

creating_table(gal_number,cloud_list,df_basic,date)
