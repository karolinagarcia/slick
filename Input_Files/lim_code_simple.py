date = '20220512'
#basicfilename = '/home/karolina.garcia/my_scripts/SLURM/BasicCharacteristics_z=1.007.csv'
#basicfilename = '/home/karolina.garcia/my_scripts/SLURM/BasicCharacteristics_z=0.0.csv'
basicfilename = '/home/karolina.garcia/my_scripts/SLURM/BasicCharacteristics_z=0.993.csv'

import pandas as pd
import argparse
from limfunctions import densityProfile, submm_luminosity, creating_table

df_basic = pd.read_csv(basicfilename)

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-ci","--cloudinfo", type=int, nargs="+")

args = parser.parse_args()

gal_number = args.cloudinfo[0]
cloud_list = args.cloudinfo[1:]

df = pd.DataFrame({'Galaxy_ID':[], 'Cloud_ID':[], 'Mcloud':[], 'Rcloud':[], 'Metallicity':[], 'RadField':[], 'redshift':[], 'H2_lcii':[], 'CO10':[], 'CO21':[], 'CO32':[], 'CO43':[], 'CO54':[], 'CO10_intTB':[], 'CO21_intTB':[], 'CO32_intTB':[], 'CO43_intTB':[], 'CO54_intTB':[], 'CI10':[], 'CI21':[], 'CO65':[], 'CO76':[], 'CO87':[], 'CO98':[], 'CO65_intTB':[], 'CO76_intTB':[], 'CO87_intTB':[], 'CO98_intTB':[], 'OI1':[], 'OI2':[], 'OI3':[], 'fH2':[]})

df.to_csv('/orange/narayanan/karolina.garcia/SLURM_'+date+'/lim_df_'+date+'.csv', index = False)

creating_table(gal_number,cloud_list,df_basic,date)