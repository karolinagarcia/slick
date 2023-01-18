# how it works

<br>

1) It takes a [SIMBA](http://simba.roe.ac.uk/) snapshot and its galaxy catalog (given by *[caesar](https://caesar.readthedocs.io/en/latest/)*), and writes the relevant information from all its particles in a pandas dataframe (which will be saved as "Basic_Characteristics*")

<br>

2) It creates a file with lists of cloud IDs (which will be saved as "Clouds_per_Core*"). Each line represents the IDs that will be passed to each core.

<br>

3) It runs radiative transfer equations on all the particles (or on a random set of those depending on the given parameters), and outputs a table with their CO, [CI] and [CII] luminosities.