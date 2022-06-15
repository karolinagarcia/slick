import numpy as np
import caesar
import yt
import pandas as pd
import yt.units as u
from tqdm import tqdm

def create_basic_table(config):
    yt_snap = yt.load(config["ytfilename"])
    yt_data = yt_snap.all_data()
    obj = caesar.load(config["caesarfilename"])

    clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun'),gal.masses['H2'].in_units('Msun'),gal.radii['baryon_half_mass'].in_units('pc'),gal.metallicities['mass_weighted']) for gal in obj.galaxies]

    kB = 1.3807 * 10**(-16) * u.cm * u.cm * u.g / (u.s * u.s * u.K)
    mH = 1.6733 * 1e-24 * u.g

    df = pd.DataFrame({'g_Index':[], 'c_Index':[], 'c_Mass':[], 'c_Radius':[], 'c_nDensity':[],
                       'c_Temperature':[], 'c_Pressure':[], 'c_Metallicity':[], 'g_SFR':[], 'g_Redshift':[],
                       'g_Mass_Gas':[], 'g_Mass_H2':[], 'g_Radius':[], 'g_Metallicity':[]})

    #for g in tqdm(np.arange(5)):
    for g in tqdm(np.arange(len(clouds_in_each_galaxy))):

        clouds_in_this_galaxy = list(clouds_in_each_galaxy[g][1])

        gal_index = [int(clouds_in_each_galaxy[g][0])]*len(clouds_in_this_galaxy)
        sfr_gal = [np.sum(yt_data['PartType0','StarFormationRate'][clouds_in_this_galaxy].value)]*len(clouds_in_this_galaxy)

        Mgas_gal = [clouds_in_each_galaxy[g][2]]*len(clouds_in_this_galaxy)
        MH2_gal = [clouds_in_each_galaxy[g][3]]*len(clouds_in_this_galaxy)
        R_gal = [clouds_in_each_galaxy[g][4]]*len(clouds_in_this_galaxy)
        Metal_gal = [clouds_in_each_galaxy[g][5]]*len(clouds_in_this_galaxy)

        Mcloud = yt_data['PartType0','Masses'][clouds_in_this_galaxy].in_units('Msun')
        n_density = yt_data['PartType0', 'Density'][clouds_in_this_galaxy].in_units('g/cm**3')/mH
        temp = yt_data['PartType0', 'Temperature'][clouds_in_this_galaxy]
        P = n_density*kB*temp
        Rcloud =  (P/(kB*1.e4))**(-0.25) * (Mcloud/(290*u.Msun))**0.5 * u.pc
        Metallicity = yt_data['PartType0','Metallicity_00'][clouds_in_this_galaxy]/.0196
        redshift = [yt_snap.parameters['Redshift']]*len(clouds_in_this_galaxy)

        df_aux = pd.DataFrame({'g_Index':gal_index, 'c_Index':clouds_in_this_galaxy, 'c_Mass':Mcloud, 'c_Radius':Rcloud,
                           'c_nDensity':n_density, 'c_Temperature':temp, 'c_Pressure':P, 'c_Metallicity':Metallicity,
                           'g_SFR':sfr_gal, 'g_Redshift':redshift, 'g_Mass_Gas':Mgas_gal, 'g_Mass_H2':MH2_gal, 'g_Radius':R_gal, 'g_Metallicity':Metal_gal})
        df = df.append(df_aux)

    df = df.reset_index(drop=True)
    df['g_Index'] = df['g_Index'].astype('int')
    df['c_Index'] = df['c_Index'].astype('int')

    basic_filename = f'Basic_Characteristics_m{config["boxsize"]}_z={round(yt_snap.parameters["Redshift"],3)}.csv'
    df.to_csv(basic_filename, index = False)
    
    #return basic_filename
