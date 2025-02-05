import numpy as np
import caesar
import yt
import pandas as pd
import yt.units as u
from tqdm import tqdm
from scipy.spatial import KDTree
import ast


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

    # basic_filename = f'/blue/narayanan/karolina.garcia/github/slick/basic_tables/Basic_Characteristics_m{config["boxsize"]}_z={round(yt_snap.parameters["Redshift"],3)}_gal5.csv'
    basic_filename = f"{config['basictable_dir']}/Basic_Characteristics_m{config['boxsize']}_z={round(yt_snap.parameters['Redshift'], 3)}_gal5.csv"
    # basic_filename = f'{config["basictable_dir"]}/Basic_Characteristics_m{config["boxsize"]}_z={round(yt_snap.parameters["Redshift"],3)}.csv'
    # basic_filename = f'{config["basictable_dir"]}/Basic_Characteristics_m{config["boxsize"]}_run={str(config["caesarfilename"])[-41:-39]}.csv'

    # removing duplicates, keeping only the first, just to make sure
    df = df.drop_duplicates(subset="c_Index")

    df.to_csv(basic_filename, index=False)

    return df
