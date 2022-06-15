import pandas as pd
import argparse
import glob

from limfunctions import creating_table
from read_config import parse_parameters

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-p", "--parameters", type=str)
parser.add_argument("--cloudinfofile", type=str)
parser.add_argument("--cloudinfoline", type=int)

args = parser.parse_args()

cloud_info_file = args.cloudinfofile
print(cloud_info_file)
with open(cloud_info_file, "r") as f:
    [gal_number, *cloud_list] = [int(x) for x in next(x for i, x in enumerate(f) if i == args.cloudinfoline).split()]

config = parse_parameters(args.parameters)

# Creates final luminosity table
df_basic = pd.read_csv(glob.glob('Output_Files/Basic*')[0])

df = pd.DataFrame({'Galaxy_ID':[], 'Cloud_ID':[], 'Mcloud':[], 'Rcloud':[], 'Metallicity':[], 'RadField':[], 'redshift':[], 'H2_lcii':[], 'CO10':[], 'CO21':[], 'CO32':[], 'CO43':[], 'CO54':[], 'CO10_intTB':[], 'CO21_intTB':[], 'CO32_intTB':[], 'CO43_intTB':[], 'CO54_intTB':[], 'CI10':[], 'CI21':[], 'CO65':[], 'CO76':[], 'CO87':[], 'CO98':[], 'CO65_intTB':[], 'CO76_intTB':[], 'CO87_intTB':[], 'CO98_intTB':[], 'OI1':[], 'OI2':[], 'OI3':[], 'fH2':[]})

# df.to_csv(f'Output_Files/lim_df_{config["date"]}.csv', index = False)

creating_table(gal_number, cloud_list, df_basic, config["date"])
