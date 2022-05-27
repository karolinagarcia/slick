# SLICK -- Suite of Line Intensity Calculations (by Karol ;)
#### SLICK is a tool that calculates many different line intensities from large sets of clouds/galaxies.
#### It was designed to work for SIMBA simulation datasets, but it will be expanded and generalized in the near future.

## How it works
#### 1) It takes a SIMBA snapshot and its caesar catalog, and write the relevant information from all its particles in a pandas dataframe (which will be saved as "Basic_Characteristics*")
#### 2) It creates a file with lists of cloud IDs (which will be saved as "Clouds_per_Core*"). Each line represents the IDs that will be passed to each core.
#### 3) It runs despotic on all the particles (or on a random set of those depending on the given parameters), and outputs a table with their CO, C+ and CII luminosities.

## How to use
#### Go to slick.py and set the parameters of interest
#### The parameters are:
#### date: date that you're running the code so that the folder and the files can be named accordingly (ex: '20220527')
#### boxsize: in the case of SIMBA, 400 or 1000 Mpc/h (ex: 400)
#### ytfilename: the path and name to the yt SIMBA snapshot that you'd like to use (ex: '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5')
#### caesarfilename: the path and name to the caesar catalog related to the yt file above (ex: '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5')
#### mode = 'total' to run code in all the clouds, 'randomize' to take a random sample (ex: 'randomize')
#### --> In the case of 'randomize', one should also set the number of galaxies they would like in this sample (ex: n_galaxies_sample = 30), and a mass range (ex: min_mass = 10^11 & max_mass = 10^12)

#### I'm working on creating a parameter file separated from the main code.
