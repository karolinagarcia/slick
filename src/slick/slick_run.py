import pandas as pd
import argparse
import glob

from .limfunctions import creating_table
from .read_config import parse_parameters
from .create_cloudspercore_list import create_cloudspercore_list
from .create_jobscript import create_jobscript
from .create_final_galaxies_table import create_gal_table


def run(cloud_info_file: str, cloudinfoline: int, parameters: str):
    with open(cloud_info_file, "r") as f:
        [*cloud_list] = [
            int(x)
            for x in next(x for i, x in enumerate(f) if i == cloudinfoline - 1).split()
        ]

    config = parse_parameters(parameters)

    #df_basic = pd.read_csv(glob.glob(f"{config['basictable_dir']}/Basic*")[0])

    if config["sim"] == 'SIMBA':
        df_basic = pd.read_csv(glob.glob(f"{config["basictable_dir"]}/Basic_Characteristics_SIMBA*{config["output_name"]}.csv")[0])
    else: df_basic = pd.read_csv(f"{config["basictable_dir"]}/Basic_Characteristics_{config["sim_code"]}_snap{config["snap_number"]}_{config["output_name"]}.csv")

    df_lim = creating_table(
        cloud_list, df_basic, config["output_dir"], int(config["n_zones"], config)
    )

    # Running remaining clouds in case any core time was not enough
    set1 = set(df_basic["c_Index"])
    set2 = set(df_lim["Cloud_ID"].unique())
    unique_to_list1 = set1 - set2
    if len(unique_to_list1) > 0:
        n_clouds_per_core = config["n_clouds_per_core"] / 10
        param_filename, max_lines = create_cloudspercore_list(config, n_clouds_per_core)
        create_jobscript(
            param_filename,
            max_lines,
            int(config["max_cores"]),
            parameters,
            config["sbatch"],
            config["module"],
        )
        with open(cloud_info_file, "r") as f:
            [*cloud_list] = [
                int(x)
                for x in next(
                    x for i, x in enumerate(f) if i == cloudinfoline - 1
                ).split()
            ]
        creating_table(
            cloud_list, df_basic, config["output_dir"], int(config["n_zones"])
        )

    # creating a galaxies table (combines the cloud luminosities of each galaxy)
    create_gal_table(config)

    print("Slick run completed successfully")
