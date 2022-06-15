import pandas as pd
import argparse
import glob

from limfunctions import creating_table
from read_config import parse_parameters

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-ci","--cloudinfo", type=int, nargs="+")
parser.add_argument("-p", "--parameters", type=str, nargs=1)

args = parser.parse_args()

gal_number = args.cloudinfo[0]
cloud_list = args.cloudinfo[1:]

config = parse_parameters(args.parameters)

# Creates final luminosity table
df_basic = pd.read_csv(glob.glob('Basic*')[0])

df = pd.DataFrame({'Galaxy_ID':[], 'Cloud_ID':[], 'Mcloud':[], 'Rcloud':[], 'Metallicity':[], 'RadField':[], 'redshift':[], 'H2_lcii':[], 'CO10':[], 'CO21':[], 'CO32':[], 'CO43':[], 'CO54':[], 'CO10_intTB':[], 'CO21_intTB':[], 'CO32_intTB':[], 'CO43_intTB':[], 'CO54_intTB':[], 'CI10':[], 'CI21':[], 'CO65':[], 'CO76':[], 'CO87':[], 'CO98':[], 'CO65_intTB':[], 'CO76_intTB':[], 'CO87_intTB':[], 'CO98_intTB':[], 'OI1':[], 'OI2':[], 'OI3':[], 'fH2':[]})

df.to_csv(f'lim_df_{config["date"]}.csv', index = False)

creating_table(gal_number,cloud_list,df_basic)
