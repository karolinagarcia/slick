import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.notebook import tqdm
from scipy.stats import binned_statistic
from scipy.integrate import quad
from colossus.cosmology import cosmology

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing


def table_for_plots(df, df_basic):
    index_99 = df[df["CO10"] == -99.0].index.tolist()
    df_no99 = df.drop(index_99)
    df_no99["SurfaceDens"] = df_no99["Mcloud"] / (np.pi * df_no99["Rcloud"] ** 2)
    df_no99["CO10_intTB"] = df_no99["CO10_intTB"] * (np.pi * df_no99["Rcloud"] ** 2)
    df_no99

    df_basic = df_basic.rename(columns={"g_Index": "Galaxy_ID"})

    j = {
        "g_SFR": "first",
        "g_Redshift": "first",
        "g_Mass_Gas": "first",
        "g_Mass_H2": "first",
        "g_Radius": "first",
        "g_Metallicity": "first",
    }
    df_gals_aux1 = df_basic.groupby(["Galaxy_ID"], as_index=False).agg(j)

    f = {
        "H2_lcii": "sum",
        "CO10": "sum",
        "CO21": "sum",
        "CO32": "sum",
        "CO43": "sum",
        "CO54": "sum",
        "CO10_intTB": "sum",
        "CO21_intTB": "mean",
        "CO32_intTB": "mean",
        "CO43_intTB": "mean",
        "CO54_intTB": "mean",
        "CI10": "sum",
        "CI21": "sum",
        "CO65": "sum",
        "CO76": "sum",
        "CO87": "sum",
        "CO98": "sum",
        "CO65_intTB": "mean",
        "CO76_intTB": "mean",
        "CO87_intTB": "mean",
        "CO98_intTB": "mean",
        "OI1": "sum",
        "OI2": "sum",
        "OI3": "sum",
        "fH2": "sum",
    }
    df_gals_aux2 = df_no99.groupby(["Galaxy_ID"], as_index=False).agg(f)

    df_gals = pd.merge(df_gals_aux1, df_gals_aux2, on="Galaxy_ID")

    return df_gals


def ToTheEnd(df, columns):
    for c in columns:
        Target_data = df[c]
        df = df.drop([c], axis=1)
        df[c] = Target_data
    return df


def train_model(df, df_basic, line):
    data = table_for_plots(df, df_basic)
    if type(df_basic["g_Mass_Gas"].iloc[0]) != np.float64:
        df_basic["g_Mass_Gas"] = (
            df_basic["g_Mass_Gas"].map(lambda x: x.rstrip(" Msun")).astype(float)
        )
        df_basic["g_Mass_H2"] = (
            df_basic["g_Mass_H2"].map(lambda x: x.rstrip(" Msun")).astype(float)
        )
        df_basic["g_Radius"] = (
            df_basic["g_Radius"].map(lambda x: x.rstrip(" pc")).astype(float)
        )
        df_basic["g_Metallicity"] = (
            df_basic["g_Metallicity"]
            .map(lambda x: x.rstrip(" dimensionless"))
            .astype(float)
        )

    input_cols = ["g_SFR", "g_Mass_Gas", "g_Radius", "g_Metallicity"]
    output_cols = [
        "H2_lcii",
        "CO10",
        "CO21",
        "CO32",
        "CO43",
        "CO54",
        "CO10_intTB",
        "CI10",
        "CI21",
        "fH2",
    ]

    # output_cols = ['H2_lcii','CO10','CO21','CO32','CO43','CO54','CO10_intTB','CO21_intTB','CO32_intTB','CO43_intTB','CO54_intTB','CI10','CI21','CO65',
    #            'CO76','CO87','CO98','CO65_intTB','CO76_intTB','CO87_intTB','CO98_intTB','OI1','OI2', 'OI3', 'fH2']
    data = ToTheEnd(data, output_cols)

    X = data.loc[:, data.columns.isin(input_cols)]
    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = scaler.transform(X)

    forest = RandomForestRegressor(1000)

    y = data[line]
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.20, random_state=44, shuffle=True
    )
    forest.fit(X_train, y_train)
    r2 = forest.score(X_test, y_test)
    # predictions = forest.predict(X_test)
    # print('\n',line,': ', round(r2*100, 0), '%')

    return forest


def predict(df_basic, line, model):
    df_basic_final = pd.DataFrame()
    df_basic = df_basic.rename(columns={"g_Index": "Galaxy_ID"})
    j = {
        "g_SFR": "first",
        "g_Redshift": "first",
        "g_Mass_Gas": "first",
        "g_Mass_H2": "first",
        "g_Radius": "first",
        "g_Metallicity": "first",
    }
    df_basic = df_basic.groupby(["Galaxy_ID"], as_index=False).agg(j)
    df_basic_final = pd.concat([df_basic_final, df_basic])

    # removing units from columns
    if type(df_basic_final["g_Mass_Gas"].iloc[0]) != np.float64:
        df_basic_final["g_Mass_Gas"] = (
            df_basic_final["g_Mass_Gas"].map(lambda x: x.rstrip(" Msun")).astype(float)
        )
        df_basic_final["g_Mass_H2"] = (
            df_basic_final["g_Mass_H2"].map(lambda x: x.rstrip(" Msun")).astype(float)
        )
        df_basic_final["g_Radius"] = (
            df_basic_final["g_Radius"].map(lambda x: x.rstrip(" pc")).astype(float)
        )
        df_basic_final["g_Metallicity"] = (
            df_basic_final["g_Metallicity"]
            .map(lambda x: x.rstrip(" dimensionless"))
            .astype(float)
        )

    input_cols = ["g_SFR", "g_Mass_Gas", "g_Radius", "g_Metallicity"]
    X = df_basic_final.loc[:, df_basic_final.columns.isin(input_cols)]
    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = scaler.transform(X)

    # Use the forest's predict method on the "real" data
    predictions_real = model.predict(X_scaled)

    df_basic_final[line + "_pred"] = predictions_real

    # df_basic_final.to_csv('lim_df_z='+str(z)+'_'+line+'_predictedfromold.csv', index = False)
    # print("Table saved on "+'lim_df_z='+str(z)+'_'+line+'_predictedfromold.csv')

    # return predictions_real
    return df_basic_final, predictions_real
