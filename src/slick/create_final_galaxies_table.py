import pandas as pd

# df_basic = pd.read_csv('../output_100Mpc_z2_16zones_selectedgalaxies_3/Basic_Characteristics_m100_z=2.025_gal5.csv')
df_basic = pd.read_csv(
    "../output_z1_16zones_total/Basic_Characteristics_m25_run=51.csv"
)
df_basic = df_basic.rename(columns={"g_Index": "Galaxy_ID"})

# df_lim = pd.read_csv('../output_100Mpc_z2_16zones_selectedgalaxies_3/lim_df.csv')
df_lim = pd.read_csv("../output_z1_16zones_total/lim_df_completewithML.csv")
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

# df_gals.to_csv('../output_100Mpc_z2_16zones_selectedgalaxies_3/gal_df.csv',index=False)
df_gals.to_csv("../output_z1_16zones_total/gal_df.csv", index=False)
