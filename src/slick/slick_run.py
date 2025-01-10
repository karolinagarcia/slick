import pandas as pd
import argparse
import glob

from limfunctions import creating_table
from read_config import parse_parameters
from create_cloudspercore_list import create_cloudspercore_list
from create_jobscript import create_jobscript

parser = argparse.ArgumentParser(prog='lim_code_simple')
parser.add_argument("-p", "--parameters", type=str)
parser.add_argument("--cloudinfofile", type=str)
parser.add_argument("--cloudinfoline", type=int)

args = parser.parse_args()

cloud_info_file = args.cloudinfofile
with open(cloud_info_file, "r") as f:
    [*cloud_list] = [int(x) for x in next(x for i, x in enumerate(f) if i == args.cloudinfoline-1).split()]

config = parse_parameters(args.parameters)

df_basic = pd.read_csv(glob.glob(f'{config["basictable_dir"]}/Basic*')[0])

df_lim = creating_table(cloud_list, df_basic, config["output_dir"], int(config["n_zones"]))

# Running remaining clouds in case any core time was not enough
set1 = set(df_basic['c_Index'])
set2 = set(df_lim['Cloud_ID'].unique())
unique_to_list1 = set1 - set2
if unique_to_list1>0:
    n_clouds_per_core = config["n_clouds_per_core"]/10
    param_filename, max_lines = create_cloudspercore_list(config,n_clouds_per_core)
    create_jobscript(param_filename, max_lines, int(config["max_cores"]), args.parameters, config["sbatch"])
    with open(cloud_info_file, "r") as f:
        [*cloud_list] = [int(x) for x in next(x for i, x in enumerate(f) if i == args.cloudinfoline-1).split()]
    creating_table(cloud_list, df_basic, config["output_dir"], int(config["n_zones"]))

print("Slick run completed successfully")