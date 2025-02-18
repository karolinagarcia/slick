import pandas as pd
import glob

def create_gal_table(config):
    if config["sim"] == 'SIMBA':
        df_basic = pd.read_csv(glob.glob(f"{config["basictable_dir"]}/Basic_Characteristics_SIMBA*{config["output_name"]}.csv")[0])
        df_lim = pd.read_csv(f"{config["output_dir"]}/Luminosity_Clouds_SIMBA{config["boxsize"]}_{config["output_name"]}.csv")
    else:
        df_basic = pd.read_csv(f"{config["basictable_dir"]}/Basic_Characteristics_{config["sim_code"]}_snap{config["snap_number"]}_{config["output_name"]}.csv")
        df_lim = pd.read_csv(f"{config["output_dir"]}/Luminosity_Clouds_{config["sim_code"]}_snap{config["snap_number"]}_{config["output_name"]}.csv")
    
    df_basic = df_basic.rename(columns={"g_Index": "Galaxy_ID"})

    # df_lim = df_lim[df_lim['NZONES']==32.0]
    df_lim["H2_Mass"] = df_lim["fH2"] * df_lim["Mcloud"]

    j = {
        "g_Redshift": "first",
        "g_Mass_Gas": "first",
        "g_Radius": "first",
        "g_Metallicity": "first",
        "g_SFR": "first",
    }
    df_gals_aux1 = df_basic.groupby(["Galaxy_ID"], as_index=False).agg(j)

    f = {
        "H2_Mass": "sum",
        "CO10": "sum",
        "H2_lcii": "sum",
        "CO10": "sum",
        "CO21": "sum",
        "CO32": "sum",
        "CO43": "sum",
        "CO54": "sum",
        "CI10": "sum",
        "CI21": "sum",
        "CO65": "sum",
        "CO76": "sum",
        "CO87": "sum",
        "CO98": "sum",
        "OI1": "sum",
        "OI2": "sum",
        "OI3": "sum",
        "fH2": "sum",
    }
    df_gals_aux2 = df_lim.groupby(["Galaxy_ID"], as_index=False).agg(f)

    df_gals = pd.merge(df_gals_aux1, df_gals_aux2, on="Galaxy_ID")

    if config["sim"] == 'SIMBA':
        # df_gals.to_csv('../output_100Mpc_z2_16zones_selectedgalaxies_3/gal_df.csv',index=False)
        df_gals.to_csv(f"{config["output_dir"]}/Luminosity_Galaxies_SIMBA{config["boxsize"]}_{config["output_name"]}.csv", index=False)
    else:
        df_gals.to_csv(f"{config["output_dir"]}/Luminosity_Galaxies_{config["sim_code"]}_snap{config["snap_number"]}_{config["output_name"]}.csv", index=False)
