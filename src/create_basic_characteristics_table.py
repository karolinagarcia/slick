import numpy as np
import caesar
import yt
import pandas as pd
import yt.units as u
from tqdm import tqdm
from scipy.spatial import KDTree

def create_basic_table(config):

    yt_snap = yt.load(config["ytfilename"])
    yt_data = yt_snap.all_data()
    obj = caesar.load(config["caesarfilename"])

    clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun').value,gal.masses['H2'].in_units('Msun').value,gal.radii['baryon_half_mass'].in_units('pc').value,gal.metallicities['mass_weighted'].value) for gal in obj.galaxies]

    kB = 1.3807 * 10**(-16) * u.cm * u.cm * u.g / (u.s * u.s * u.K)
    mH = 1.6733 * 1e-24 * u.g
    K_abs = 1.07800e5 * u.cm * u.cm / u.g   #median value of the mass attenuation coefficient (absorption cross section per g of dust) in cm^2/g within the Habing limit of 91.2 nm to 111.0 nm (from the Draine table)

    df = pd.DataFrame({'g_Index':[], 'c_Index':[], 'c_Mass':[], 'c_Radius':[], 'c_nDensity':[],
                       'c_Temperature':[], 'c_Pressure':[], 'c_Metallicity':[], 'c_SFR':[], 'c_RadField':[], 'c_DMR':[], 'g_SFR':[], 'g_Redshift':[],
                       'g_Mass_Gas':[], 'g_Mass_H2':[], 'g_Radius':[], 'g_Metallicity':[]})

    for cloud in tqdm(clouds_in_each_galaxy):

        clouds_in_this_galaxy = list(cloud[1])     #gives the list of indices of gas particles in the galaxy indexed by g. len(clouds_in_this_galaxy) thus gives the number of clouds in the specified galaxy

        gal_index = [int(cloud[0])]*len(clouds_in_this_galaxy)
        sfr_gal = [np.sum(yt_data['PartType0','StarFormationRate'][clouds_in_this_galaxy].value)]*len(clouds_in_this_galaxy)

        Mgas_gal = [cloud[2]]*len(clouds_in_this_galaxy)    #mass of gas in the galaxy indexed by g duplicated to a list of len(clouds_in_this_galaxy) elements
        MH2_gal = [cloud[3]]*len(clouds_in_this_galaxy)
        R_gal = [cloud[4]]*len(clouds_in_this_galaxy)
        Metal_gal = [cloud[5]]*len(clouds_in_this_galaxy)

        Mcloud = yt_data['PartType0','Masses'][clouds_in_this_galaxy].in_units('Msun')    #making an array of each cloud's masses indexed according to glist
        n_density = yt_data['PartType0', 'Density'][clouds_in_this_galaxy].in_units('g/cm**3')/mH    #making an array of each cloud's number densities
        temp = yt_data['PartType0', 'Temperature'][clouds_in_this_galaxy]    #making an array of each cloud's temperatures
        P = n_density*kB*temp
        Rcloud =  (P/(kB*1.e4*u.K/(u.cm*u.cm*u.cm)))**(-0.25) * (Mcloud/(290*u.Msun))**0.5 * u.pc    #making an array of each cloud's radius
        Metallicity = yt_data['PartType0','Metallicity_00'][clouds_in_this_galaxy]/.0196    #making an array of each cloud's metallicity normalized to the solar neighbourhood (MW) (based on solar wind metallicity)
        redshift = [yt_snap.parameters['Redshift']]*len(clouds_in_this_galaxy)

        #identifying all the neighbours first

        coords_gas = yt_data["PartType0","Coordinates"][clouds_in_this_galaxy].in_units("kpc")
        
        gas_dust = yt_snap.arr(yt_data['PartType0','Dust_Masses'][clouds_in_this_galaxy],"code_mass").in_units('Msun')
    
        gas_kd_tree = KDTree(coords_gas.value)
        gas_tree_dist, gas_indexes = gas_kd_tree.query(coords_gas.value, k=64)    #finding 64 nearest neighbour gas particles, dist in units of kpc

        sfr_gas = yt_data["PartType0","StarFormationRate"][clouds_in_this_galaxy].in_units("Msun/yr")

        RadField_list = []
        Metallicity_val = cloud[5]
        if Metallicity_val==0:
            DMR = 0
        else:
            DGR = (10**(2.445*np.log10(Metallicity_val/0.0134)-2.029))    #equation 9, Qi Li+2019
            DMR = DGR/(Metallicity_val*0.44)   #corrected DMR = dust mass/metal mass = dust mass/(metallicity*total mass). This is normalized to the MW DMR value of 0.44        
        DMR_list = [DMR]*len(clouds_in_this_galaxy)

        
        for i in range(0,len(clouds_in_this_galaxy)):
            try:
                cross_section_radius = np.sqrt(((coords_gas[i][0]-coords_gas[gas_indexes[i][-1]][0]).value)**2+((coords_gas[i][1]-coords_gas[gas_indexes[i][-1]][1]).value)**2+((coords_gas[i][2]-coords_gas[gas_indexes[i][-1]][2]).value)**2) * u.kpc    #distance to 64th particle in units of kpc
                sfr_surface_density = np.sum(sfr_gas[gas_indexes[i]].value)/(np.pi*(cross_section_radius.value)**2)    #this is in units of Msun/yr/kpc**2
                solar_sfr_surface_density = 790e-6    #in units of Msun/yr/kpc**2 ; value from Bonatto & Bica, 2011
            
                optical_depth = K_abs * (cross_section_radius.in_units("cm")) * (np.sum(gas_dust[gas_indexes[i]])/(4/3*np.pi*(cross_section_radius)**3)).in_units("g/cm**3")
            
                if optical_depth==0:
                    beta_UV = 1
                else:
                    beta_UV = (1-np.exp(-optical_depth.value))/(optical_depth.value)    #equation from Claudia Lagos 2014
            
                solar_gas_surface_density = 10 * u.Msun / (u.pc * u.pc)    #Chang et al 2002 (Scoville & Sanders 1987)
                solar_optical_depth = K_abs * solar_gas_surface_density.in_units("g/cm**2")/1.653e2    #gas to dust ratio used by the Draine table used here 
                solar_beta = (1-np.exp(-solar_optical_depth.value))/(solar_optical_depth.value)
            
                RadField_val = (sfr_surface_density/solar_sfr_surface_density)*(beta_UV/solar_beta)
                RadField_list.append(RadField_val)
            except:
                RadField_list.append(0)

        df_aux = pd.DataFrame({'g_Index':gal_index, 'c_Index':clouds_in_this_galaxy, 'c_Mass':Mcloud, 'c_Radius':Rcloud, 'c_nDensity':n_density, 'c_Temperature':temp, 'c_Pressure':P, 'c_Metallicity':Metallicity, 'c_SFR':sfr_gas, 'c_RadField':RadField_list, 'c_DMR':DMR_list, 'g_SFR':sfr_gal, 'g_Redshift':redshift, 'g_Mass_Gas':Mgas_gal, 'g_Mass_H2':MH2_gal, 'g_Radius':R_gal, 'g_Metallicity':Metal_gal})
        df = pd.concat([df, df_aux])

    df = df.reset_index(drop=True)
    df['g_Index'] = df['g_Index'].astype('int')
    df['c_Index'] = df['c_Index'].astype('int')

    basic_filename = f'../{config["basictable_dir"]}/Basic_Characteristics_m{config["boxsize"]}_z={round(yt_snap.parameters["Redshift"],3)}.csv'
    df.to_csv(basic_filename, index = False)
