import numpy as np
import caesar
import yt
import pandas as pd
import yt.units as u
from tqdm import tqdm
from scipy.spatial import KDTree
import astropy.constants as ac
import h5py
from colossus.cosmology import cosmology
import requests


### SIMBA

def create_basic_table(config):
    yt_snap = yt.load(config["ytfilename"])
    yt_data = yt_snap.all_data()
    obj = caesar.load(config["caesarfilename"])

    clouds_in_each_galaxy = [
        (
            gal.GroupID,
            gal.glist,
            gal.masses["gas"].in_units("Msun").value,
            gal.masses["H2"].in_units("Msun").value,
            gal.radii["gas_half_mass"].in_units("pc").value,
            gal.metallicities["mass_weighted"].value,
            gal.pos.in_units("kpc").value,
        )
        for gal in obj.galaxies
    ]

    # print(config['agn'])
    """
    if config['agn'] == 'True':
        clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun').value,gal.masses['H2'].in_units('Msun').value,gal.radii['gas_half_mass'].in_units('pc').value,gal.metallicities['mass_weighted'].value,gal.pos.in_units('kpc').value,gal.masses['bh'].in_units('Msun').value,gal.bhmdot.in_units('g/s').value,gal.bh_fedd.value) for gal in obj.galaxies]
    # this should be named galaxies

    if config["mode"] == 'single':
        clouds_in_each_galaxy = [clouds_in_each_galaxy[i] for i in ast.literal_eval(config['gal_ids'])]
    """
    c = 2.99792458 * 1e10 * u.cm / u.s
    L_sun = 3.826 * 10**33 * u.erg / u.s
    kB = 1.3807 * 10 ** (-16) * u.cm * u.cm * u.g / (u.s * u.s * u.K)
    mH = 1.6733 * 1e-24 * u.g
    K_abs = (
        1.07800e5 * u.cm * u.cm / u.g
    )  # median value of the mass attenuation coefficient (absorption cross section per g of dust) in cm^2/g within the Habing limit of 91.2 nm to 111.0 nm (from the Draine table)
    """
    if config['agn'] == 'True':
        df = pd.DataFrame({'g_Index':[], 'c_Index':[], 'c_Mass':[], 'c_Radius':[], 'c_nDensity':[],
                       'c_Temperature':[], 'c_Pressure':[], 'c_Metallicity':[], 'c_SFR':[], 'c_RadField':[], 'c_CR':[] ,'c_DMR':[], 'g_SFR':[], 'g_Redshift':[],
                       'g_Mass_Gas':[], 'g_Mass_H2':[], 'g_Radius':[], 'g_Metallicity':[],
                       'BHmass_gal':[], 'BHfedd_gal':[], 'L_edd':[], 'BH_bol_lum':[],'bol_corr_UV':[],
                       'BH_UV_lum':[],'BH_bol_lum_2':[],'BH_UV_lum_2':[],'BH_UV_flux':[]})

    else:
        df = pd.DataFrame({'g_Index':[], 'c_Index':[], 'c_Mass':[], 'c_Radius':[], 'c_nDensity':[],
                       'c_Temperature':[], 'c_Pressure':[], 'c_Metallicity':[], 'c_SFR':[], 'c_RadField':[], 'c_CR':[] ,'c_DMR':[], 'g_SFR':[], 'g_Redshift':[],
                       'g_Mass_Gas':[], 'g_Mass_H2':[], 'g_Radius':[], 'g_Metallicity':[]})

    """
    df = pd.DataFrame(
        {
            "g_Index": [],
            "c_Index": [],
            "c_Mass": [],
            "c_Radius": [],
            "c_nDensity": [],
            "c_Temperature": [],
            "c_Pressure": [],
            "c_Metallicity": [],
            "c_SFR": [],
            "c_RadField": [],
            "c_CR": [],
            "c_DMR": [],
            "g_SFR": [],
            "g_Redshift": [],
            "g_Mass_Gas": [],
            "g_Mass_H2": [],
            "g_Radius": [],
            "g_Metallicity": [],
        }
    )

    for cloud in tqdm(clouds_in_each_galaxy):
        # cloud should be named galaxy

        clouds_in_this_galaxy = list(
            cloud[1]
        )  # gives the list of indices of gas particles in the galaxy indexed by g. len(clouds_in_this_galaxy) thus gives the number of clouds in the specified galaxy

        gal_index = [int(cloud[0])] * len(clouds_in_this_galaxy)
        sfr_gal = [
            np.sum(
                yt_data["PartType0", "StarFormationRate"][clouds_in_this_galaxy].value
            )
        ] * len(clouds_in_this_galaxy)

        Mgas_gal = (
            [cloud[2]] * len(clouds_in_this_galaxy)
        )  # mass of gas in the galaxy indexed by g duplicated to a list of len(clouds_in_this_galaxy) elements
        MH2_gal = [cloud[3]] * len(clouds_in_this_galaxy)
        R_gal = [cloud[4]] * len(clouds_in_this_galaxy)
        Metal_gal = [cloud[5] / 0.0196] * len(clouds_in_this_galaxy)
        BHpos_gal = [cloud[6]] * len(clouds_in_this_galaxy)
        """
        if config['agn'] == 'True':

            BHmass_gal = cloud[7]
            BHmdot_gal = cloud[8]
            BHfedd_gal = cloud[9]

            ### CHANGE THIS
            L_edd = 3.28 * 10**4 * BHmass_gal  # in L_sun
            BH_bol_lum = BHfedd_gal * L_edd
            bol_corr_UV = 5   # from Richards et al 2006 https://iopscience.iop.org/article/10.1086/506525/pdf
            BH_UV_lum = BH_bol_lum/bol_corr_UV

            rad_eff = 0.1
            if BHfedd_gal > 0.1:
                BH_bol_lum_2 = (rad_eff/(1-rad_eff))*BHmdot_gal*c.value**2  #erg/s
            else:
                BH_bol_lum_2 = 10*BHfedd_gal*rad_eff*BHmdot_gal*c.value**2  #erg/s

            BH_UV_lum_2 = BH_bol_lum_2/bol_corr_UV

            BHmass_gal = [BHmass_gal]*len(clouds_in_this_galaxy)
            BHmdot_gal = [BHmdot_gal]*len(clouds_in_this_galaxy)
            BHfedd_gal = [BHfedd_gal]*len(clouds_in_this_galaxy)
            L_edd = [L_edd]*len(clouds_in_this_galaxy)
            BH_bol_lum = [BH_bol_lum]*len(clouds_in_this_galaxy)
            bol_corr_UV = [bol_corr_UV]*len(clouds_in_this_galaxy)
            BH_UV_lum = [BH_UV_lum]*len(clouds_in_this_galaxy)
            BH_bol_lum_2 = [BH_bol_lum_2]*len(clouds_in_this_galaxy)
            BH_UV_lum_2 = [BH_UV_lum_2]*len(clouds_in_this_galaxy)
        """

        Mcloud = yt_data["PartType0", "Masses"][clouds_in_this_galaxy].in_units(
            "Msun"
        )  # making an array of each cloud's masses indexed according to glist
        n_density = (
            yt_data["PartType0", "Density"][clouds_in_this_galaxy].in_units("g/cm**3")
            / mH
        )  # making an array of each cloud's number densities
        temp = yt_data["PartType0", "Temperature"][
            clouds_in_this_galaxy
        ]  # making an array of each cloud's temperatures
        P = n_density * kB * temp
        Rcloud = (
            (P / (kB * 1.0e4 * u.K / (u.cm * u.cm * u.cm))) ** (-0.25)
            * (Mcloud / (290 * u.Msun)) ** 0.5
            * u.pc
        )  # making an array of each cloud's radius
        Metallicity = (
            yt_data["PartType0", "Metallicity_00"][clouds_in_this_galaxy] / 0.0196
        )  # making an array of each cloud's metallicity normalized to the solar neighbourhood (MW) (based on solar wind metallicity)
        # add abundances Metallicity_1, 2, 3 etc... carbon and oxygen abundances
        #
        redshift = [yt_snap.parameters["Redshift"]] * len(clouds_in_this_galaxy)

        # identifying all the neighbours first

        coords_gas = yt_data["PartType0", "Coordinates"][
            clouds_in_this_galaxy
        ].in_units("kpc")
        # coords_BH = yt_data["PartType5","Coordinates"][clouds_in_this_galaxy].in_units("kpc")

        gas_dust = yt_snap.arr(
            yt_data["PartType0", "Dust_Masses"][clouds_in_this_galaxy], "code_mass"
        ).in_units("Msun")

        gas_kd_tree = KDTree(coords_gas.value)
        gas_tree_dist, gas_indexes = gas_kd_tree.query(
            coords_gas.value, k=64
        )  # finding 64 nearest neighbour gas particles, dist in units of kpc

        sfr_gas = yt_data["PartType0", "StarFormationRate"][
            clouds_in_this_galaxy
        ].in_units("Msun/yr")

        CR_list = []
        RadField_list = []
        UV_flux_list = []
        Metallicity_val = cloud[5]
        if Metallicity_val == 0:
            DMR = 0
        else:
            DGR = 10 ** (
                2.445 * np.log10(Metallicity_val / 0.0134) - 2.029
            )  # equation 9, Qi Li+2019
            DMR = (
                DGR / (Metallicity_val * 0.44)
            )  # corrected DMR = dust mass/metal mass = dust mass/(metallicity*total mass). This is normalized to the MW DMR value of 0.44
        DMR_list = [DMR] * len(clouds_in_this_galaxy)

        for i in range(0, len(clouds_in_this_galaxy)):
            try:
                cross_section_radius = (
                    np.sqrt(
                        ((coords_gas[i][0] - coords_gas[gas_indexes[i][-1]][0]).value)
                        ** 2
                        + ((coords_gas[i][1] - coords_gas[gas_indexes[i][-1]][1]).value)
                        ** 2
                        + ((coords_gas[i][2] - coords_gas[gas_indexes[i][-1]][2]).value)
                        ** 2
                    )
                    * u.kpc
                )  # distance to 64th particle in units of kpc
                sfr_surface_density = np.sum(sfr_gas[gas_indexes[i]].value) / (
                    np.pi * (cross_section_radius.value) ** 2
                )  # this is in units of Msun/yr/kpc**2
                solar_sfr_surface_density = 790e-6  # in units of Msun/yr/kpc**2 ; value from Bonatto & Bica, 2011

                optical_depth = (
                    K_abs
                    * (cross_section_radius.in_units("cm"))
                    * (
                        np.sum(gas_dust[gas_indexes[i]])
                        / (4 / 3 * np.pi * (cross_section_radius) ** 3)
                    ).in_units("g/cm**3")
                )

                if optical_depth == 0:
                    beta_UV = 1
                else:
                    beta_UV = (1 - np.exp(-optical_depth.value)) / (
                        optical_depth.value
                    )  # equation from Claudia Lagos 2014

                solar_gas_surface_density = (
                    10 * u.Msun / (u.pc * u.pc)
                )  # Chang et al 2002 (Scoville & Sanders 1987)
                solar_optical_depth = (
                    K_abs * solar_gas_surface_density.in_units("g/cm**2") / 1.653e2
                )  # gas to dust ratio used by the Draine table used here
                solar_beta = (1 - np.exp(-solar_optical_depth.value)) / (
                    solar_optical_depth.value
                )

                CR_val = (
                    sfr_surface_density / solar_sfr_surface_density
                )  # using the unattenuated SFR surface density for scaling the CRs
                RadField_val = (sfr_surface_density / solar_sfr_surface_density) * (
                    beta_UV / solar_beta
                )

            except:
                RadField_val = 0  # or -99

            # including BH
            BH_UV_flux = 0
            """
            if config['agn'] == 'True':
                #BHmass_treshold = 0
                #if BHmass_gal[0]>BHmass_treshold:
                part_BH_dist = np.sqrt(((coords_gas[i][0].value-BHpos_gal[0][0]))**2+((coords_gas[i][1].value-BHpos_gal[0][1]))**2+((coords_gas[i][2].value-BHpos_gal[0][2]))**2) * u.kpc    #distance to BH in units of kpc
                part_BH_dist = part_BH_dist.in_units("cm").value
                BH_UV_flux = BH_UV_lum_2[0] / (4*np.pi*part_BH_dist**2)   # in erg/s/cm**2/(sr?) --- so Habing?

                print(RadField_val)
                print(BH_UV_flux)

                RadField_val = RadField_val + BH_UV_flux
            """

            CR_list.append(CR_val)
            RadField_list.append(RadField_val)
            UV_flux_list.append(BH_UV_flux)

            # except:
            #    RadField_list.append(-99)
            #    UV_flux_list.append(-99)

        # check which galaxy this cloud belongs to
        # check if this galaxy has a SMBH
        # check distance BH to cloud
        # if whithin certain tresholds (min and max), we modify the UV rad field
        """
        if config['agn'] == 'True':
            df_aux = pd.DataFrame({'g_Index':gal_index, 'c_Index':clouds_in_this_galaxy, 'c_Mass':Mcloud,
                                'c_Radius':Rcloud, 'c_nDensity':n_density, 'c_Temperature':temp,
                                'c_Pressure':P, 'c_Metallicity':Metallicity, 'c_SFR':sfr_gas,
                                'c_RadField':RadField_list, 'c_CR':CR_list ,'c_DMR':DMR_list, 'g_SFR':sfr_gal,
                                'g_Redshift':redshift, 'g_Mass_Gas':Mgas_gal, 'g_Mass_H2':MH2_gal,
                                'g_Radius':R_gal, 'g_Metallicity':Metal_gal, 'BHmass_gal':BHmass_gal,
                                'BHfedd_gal':BHfedd_gal, 'L_edd':L_edd, 'BH_bol_lum':BH_bol_lum,
                                'BH_UV_lum':BH_UV_lum, 'BH_bol_lum_2':BH_bol_lum_2,'BH_UV_lum_2':BH_UV_lum_2,
                                'BH_UV_flux':UV_flux_list})
        else:
            df_aux = pd.DataFrame({'g_Index':gal_index, 'c_Index':clouds_in_this_galaxy, 'c_Mass':Mcloud,
                                'c_Radius':Rcloud, 'c_nDensity':n_density, 'c_Temperature':temp,
                                'c_Pressure':P, 'c_Metallicity':Metallicity, 'c_SFR':sfr_gas,
                                'c_RadField':RadField_list, 'c_CR':CR_list ,'c_DMR':DMR_list, 'g_SFR':sfr_gal,
                                'g_Redshift':redshift, 'g_Mass_Gas':Mgas_gal, 'g_Mass_H2':MH2_gal,
                                'g_Radius':R_gal, 'g_Metallicity':Metal_gal})
        """
        df_aux = pd.DataFrame(
            {
                "g_Index": gal_index,
                "c_Index": clouds_in_this_galaxy,
                "c_Mass": Mcloud,
                "c_Radius": Rcloud,
                "c_nDensity": n_density,
                "c_Temperature": temp,
                "c_Pressure": P,
                "c_Metallicity": Metallicity,
                "c_SFR": sfr_gas,
                "c_RadField": RadField_list,
                "c_CR": CR_list,
                "c_DMR": DMR_list,
                "g_SFR": sfr_gal,
                "g_Redshift": redshift,
                "g_Mass_Gas": Mgas_gal,
                "g_Mass_H2": MH2_gal,
                "g_Radius": R_gal,
                "g_Metallicity": Metal_gal,
            }
        )

        df = pd.concat([df, df_aux])

    df = df.reset_index(drop=True)
    df["g_Index"] = df["g_Index"].astype("int")
    df["c_Index"] = df["c_Index"].astype("int")

    basic_filename = f"{config['basictable_dir']}/Basic_Characteristics_SIMBA{config['boxsize']}_z={round(yt_snap.parameters['Redshift'], 3)}_{config['output_name']}.csv"
    # basic_filename = f'{config["basictable_dir"]}/Basic_Characteristics_m{config["boxsize"]}_run={str(config["caesarfilename"])[-41:-39]}.csv'

    # removing duplicates, keeping only the first, just to make sure
    df = df.drop_duplicates(subset="c_Index")

    df.to_csv(basic_filename, index=False)

    return df


### ILLUSTRIS / CAMELS

def get(path, params=None):

    # setting up base directory
    baseUrl = 'http://www.tng-project.org/api/'
    headers = {"api-key":"deda270560cd77905bc8afba45bf3143"}
    
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def create_basic_table_tng_camels(config):

    cosmo = cosmology.setCosmology('planck18')
    h = cosmo.h
    c = ac.c
    K_abs = 1.07800e5 # in cm**2/g
    K_abs_kpcMsun = K_abs * 1e-43 / 5e-34

    df = pd.DataFrame({'g_Index':[], 'c_Index': [], 'c_x': [], 'c_y': [], 'c_z': [],
            'c_Mass': [], 'c_Radius': [],
            'c_Density': [], 'c_Metallicity': [], 'c_SFR': [], 'c_RadField': [],
            'g_DMR': [], 'g_SFR': [], 'g_Redshift': [], 'g_Mass_Gas':[],
            'g_Radius':[], 'g_Metallicity':[]})
    
    if config["sim"] == 'TNG':

        # setting up what snapshot/halo we want
        # the halo we select form the code in JupyterLab
        snap_prog_url = "http://www.tng-project.org/api/"+config["sim_code"]+"/snapshots/"+config["snap_number"]+"/"
        snap_prog = get(snap_prog_url)
        gal_z = snap_prog['redshift']

        gal_ids = config["subhalo_number"].strip('[]').split(',')
        gal_ids = [int(num) for num in gal_ids]

    if config["sim"] == 'CAMELS':

        # finding the clouds of this specific galaxy
        f = h5py.File("/blue/narayanan/karolina.garcia/github/slick_CAMELS/FOF_Subfind/"+config["sim_type"]+"/"+config["sim_code"]+"/"+config["sim_code_2"]+"/fof_subhalo_tab_0"+config["snap_number"]+".hdf5", 'r')
        index       = int(config["subhalo_number"])
        len_sh      = f['Subhalo/SubhaloLen'][:]
        IDs_sh      = f['IDs/ID'][:]
        start       = np.sum(len_sh[:index])
        end         = start+len_sh[index]
        indexes_sh  = IDs_sh[start:end]
        common_indexes, indexes1, indexes2 = np.intersect1d(IDs_sh, indexes_sh, assume_unique=False, return_indices=True)


    for gal_id in gal_ids:

        sub_prog_url = snap_prog_url+"subhalos/"+str(gal_id)+"/"

        # getting halo info
        sub_prog = get(sub_prog_url)

        # printing x and y positions of this specific halo as an example
        # list of all subhalo parameters is here: https://www.tng-project.org/api/TNG50-1/snapshots/0/subhalos/2/
        sub_prog['pos_x'], sub_prog['pos_y']

        # here we get the cutout file. I don't quite understand what this file is exaclty yet,
        # but I think it access the particles info of a certain halo. we set the info we want in 'cutout_request'
        # list of parameters (physical properties of each gas particle) that could go to cutout request:
        # https://www.tng-project.org/data/docs/specifications/#parttype0
        cutout_request = {'gas':'Coordinates,Masses,GFM_Metallicity,StarFormationRate,Density,ParticleIDs,InternalEnergy,ElectronAbundance'}
        cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)

        f = h5py.File(cutout,'r')

        # extracting the properties of gas particles for the halo we specified in the beginning og this code
        part_id = f['PartType0']['ParticleIDs'][:]
        part_mass = f['PartType0']['Masses'][:] * 1e10/h * u.Msun # in M_sun
        part_dens = f['PartType0']['Density'][:] * (1e10*u.Msun/h) / (u.kpc/h)**3 # in M_sun/(kpc)**3 after all conversions
        part_rad_from_dens = ((part_mass/((4/3)*np.pi*part_dens))**(1/3)).in_units('pc') # in pc

        k_B = 1.380649e-23 * u.J / u.K
        m_p = 1.6726219e-27 * u.kg
        int_energy = f['PartType0']['InternalEnergy'][:] * (u.km / u.s)**2
        elec_abund = f['PartType0']['ElectronAbundance'][:] #non dimens.
        part_temp = ((2/3) * (int_energy/k_B) * (4/(3.28+3.04*elec_abund)) * m_p).in_units('K')
        part_ndens = (f['PartType0']['Density'][:] * (1e10*u.Msun/h) / (u.kpc/h)**3).in_units('g/cm**3')/m_p # in M_sun/(kpc)**3 after all conversions
        part_press = (part_ndens.in_units('1/cm**3')*k_B*part_temp)
        part_rad_from_temp =  (part_press/(k_B*1.e4*u.K/(u.cm*u.cm*u.cm)))**(-0.25) * (part_mass/(290*u.Msun))**0.5 * u.pc    #making an array of each cloud's radius
        #part_rad_from_sl = f['PartType0']['SubfindHsml'][:] * h #kpc

        part_metal = f['PartType0']['GFM_Metallicity'][:]/0.0127 * u.Zsun # in Z_sun after dividing
        part_dust = 0.4*part_metal*0.0127*part_mass # this is 0.4 * metal mass (which is Z absolute times particle mass)
        part_sfr = f['PartType0']['StarFormationRate'][:] * u.Msun/u.yr # in M_sun/yr
        part_coords = f['PartType0']['Coordinates'][:] / h #* u.kpc # in kpc
        part_x = f['PartType0']['Coordinates'][:,0] / h * u.kpc # in kpc
        part_y = f['PartType0']['Coordinates'][:,1] / h * u.kpc # in kpc
        part_z = f['PartType0']['Coordinates'][:,2] / h * u.kpc # in kpc

        # now extracting halo properties, and multiplying by the number of gas particles in this halo so that we populate all the lines
        gal_gas_mass = [sub_prog['mass_gas'] * 1e10/h]*len(part_id) #* u.Msun # units?
        gal_rad = [sub_prog['halfmassrad_gas'] * 1e3/h]*len(part_id) #* u.pc # units? # in pc
        gal_metal = [sub_prog['gasmetallicity']/0.0127]*len(part_id) #* u.Zsun # units? # Zsun?
        gal_sfr = [sub_prog['sfr']]*len(part_id) #* u.Msun / u.yr # in M_sun/yr

        # calculating DMR
        if gal_metal[0]==0:
            gal_DMR = [0]*len(part_id)
        else:
            gal_DGR = (10**(2.445*np.log10(gal_metal[0])-2.029)) #equation 9, Qi Li+2019   #for simba we divide by /0.0134 here to turn to solar metallicity
            gal_DMR = [gal_DGR/(gal_metal[0]*0.44)]*len(part_id) # DMR = dust mass/metal mass = dust mass/(metallicity*total mass), then normalized to the MW DMR value of 0.44

        # calculating radiation field
        gal_radfield = []
        # identifying all the neighbours first
        part_kdtree = KDTree(part_coords)
        part_kdtree_dist, part_kdtree_ids = part_kdtree.query(part_coords, k=64)    #finding 64 nearest neighbour gas particles, dist in units of kpc
        # now looping over all particles
        for i in range(0,len(part_id)):
            try:
                cross_section_radius = np.sqrt((part_coords[i][0]-part_coords[part_kdtree_ids[i][-1]][0])**2
                                                +(part_coords[i][1]-part_coords[part_kdtree_ids[i][-1]][1])**2
                                                +(part_coords[i][2]-part_coords[part_kdtree_ids[i][-1]][2])**2)
                                                # distance to 64th particle in units of kpc
                print('cross_section_radius')
                print(cross_section_radius)
                sfr_surface_density = np.sum(part_sfr[part_kdtree_ids[i]].value)/(np.pi*(cross_section_radius)**2)
                                # in Msun/yr/kpc**2
                print('sfr_surface_density')
                print(sfr_surface_density)
                solar_sfr_surface_density = 790e-6
                                # in Msun/yr/kpc**2 ; value from Bonatto & Bica, 2011
                optical_depth = K_abs_kpcMsun * cross_section_radius * np.sum(part_dust[part_kdtree_ids[i]].value)/(4/3*np.pi*(cross_section_radius)**3)
                                # non dimensional
                print('optical_depth')
                print(optical_depth)
                if optical_depth==0:
                    beta_UV = 1
                else:
                    beta_UV = (1-np.exp(-optical_depth))/(optical_depth) # non dimensional # equation from Claudia Lagos 2014
                print('beta_UV')
                print(beta_UV)
                solar_gas_surface_density_Msunpc = 10 # in Msun/pc**2 ; Chang et al 2002 (Scoville & Sanders 1987)
                solar_gas_surface_density = solar_gas_surface_density_Msunpc * 1.989e+33 / (3.086e18)**2 # in g/cm**2
                solar_optical_depth = K_abs * solar_gas_surface_density/1.653e2  # non dimensional  #gas to dust ratio used by the Draine table used here 
                solar_beta = (1-np.exp(-solar_optical_depth))/(solar_optical_depth) # non dimensional
                
                print('solar_gas_surface_density')
                print(solar_gas_surface_density)

                print('solar_beta')
                print(solar_beta)

                RadField = (sfr_surface_density/solar_sfr_surface_density)*(beta_UV/solar_beta)
                print('RadField')
                print(RadField)
                gal_radfield.append(RadField)
            except:
                gal_radfield.append(0)

        d = {'g_Index':gal_id, 'c_Index': part_id, 'c_x': part_x.value, 'c_y': part_y.value, 'c_z': part_z.value,
            'c_Mass': part_mass.value, 'c_Radius': part_rad_from_temp.value,
            'c_Density': part_dens.value, 'c_Metallicity': part_metal.value, 'c_SFR': part_sfr.value, 'c_RadField': gal_radfield,
            'g_DMR': gal_DMR, 'g_SFR': gal_sfr, 'g_Redshift': gal_z, 'g_Mass_Gas':gal_gas_mass,
            'g_Radius':gal_rad, 'g_Metallicity':gal_metal} #'c_Radius_Hsml': part_rad_hsml
        # Differences to SIMBA code:
        # I don't actually need Temperature and Pressure;
        # TNG doesn't have g_Mass_H2;
        # nDensity --> Density;
        # I added the positions
        # Two different kinds of radii

        df_aux = pd.DataFrame(data=d)
        df = pd.concat([df,df_aux],ignore_index=True)

    df['g_Index'] = df['g_Index'].astype(int)
    df['c_Index'] = df['c_Index'].astype(int)

    #basic_filename = f'{config["basictable_dir"]}/Basic_Characteristics_'+config["sim_code"]+'_snap'+config["snap_number"]+'_subhalo'+config["subhalo_number"]+'.csv'
    basic_filename = f'{config["basictable_dir"]}/Basic_Characteristics_'+config["sim_code"]+'_snap'+config["snap_number"]+'_'+config["output_name"]+'.csv'
 
    #removing duplicates, keeping only the first, just to make sure
    df = df.drop_duplicates(subset='c_Index')

    df.to_csv(basic_filename, index = False)

    return df