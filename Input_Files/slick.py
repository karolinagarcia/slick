import pandas as pd
# import seaborn as sns
# import sys
# import yt.units as u
import argparse
# from tqdm import tqdm
import random
random.seed(10)

from create_cloudspercore_list import create_param_file
from create_basic_characteristics_table import create_basic_table
# from limfunctions import densityProfile, submm_luminosity
from limfunctions import creating_table

from configparser import ConfigParser
config = ConfigParser()
config.read('parameters.ini')

date = config['date']['date']
boxsize = config['snap']['boxsize']
ytfilename = config['snap']['ytfilename']
caesarfilename = config['snap']['caesarfilename']
mode = config['sample']['mode']
n_galaxies_sample = config['sample']['n_galaxies_sample']
min_mass = config['sample']['min_mass']
max_mass = config['sample']['max_mass']

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
