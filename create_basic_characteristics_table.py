import numpy as np
import caesar
import yt
import pandas as pd
import seaborn as sns
import sys
import yt.units as u
#from tqdm.notebook import tqdm
from tqdm import tqdm
# if not notebook

#yt_snap = yt.load('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5')
#yt_snap = yt.load('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_212.hdf5')
yt_snap = yt.load('/orange/narayanan/karolina.garcia/simba/m100n1024/snaps/snap_m100n1024_105.hdf5.1')
yt_data = yt_snap.all_data()
#caesar_file = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5'
#caesar_file = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0212_z1.007.hdf5'
caesar_file = '/orange/narayanan/d.zimmerman/simba/m100n1024/caesar_cats/caesar_simba_105.hdf5'
obj = caesar.load(caesar_file)

clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun'),gal.masses['H2'].in_units('Msun'),gal.radii['baryon_half_mass'].in_units('pc'),gal.metallicities['mass_weighted']) for gal in obj.galaxies]

kB = 1.3807 * 10**(-16) * u.cm * u.cm * u.g / (u.s * u.s * u.K)
mH = 1.6733 * 1e-24 * u.g

df = pd.DataFrame({'g_Index':[], 'c_Index':[], 'c_Mass':[], 'c_Radius':[], 'c_nDensity':[],
                   'c_Temperature':[], 'c_Pressure':[], 'c_Metallicity':[], 'g_SFR':[], 'g_Redshift':[],
                   'g_Mass_Gas':[], 'g_Mass_H2':[], 'g_Radius':[], 'g_Metallicity':[]})

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

print(df)

df = df.reset_index(drop=True)
df['g_Index'] = df['g_Index'].astype('int')
df['c_Index'] = df['c_Index'].astype('int')

df.to_csv('/home/karolina.garcia/my_scripts/SLURM/BasicCharacteristics_z='+str(round(yt_snap.parameters['Redshift'],3))+'.csv', index = False)