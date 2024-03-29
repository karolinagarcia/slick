# SLICK -- Scalable Line Intensity Computation Kit

**_slick_** is a tool to calculate CO, [CI] and [CII] line intensities from extremely large datasets of clouds/galaxies in hydrodynamical simulations. It runs radiative transfer equations on each single particle of the simulation, and outputs a table with all clouds' and galaxies' luminosities for all redshifts of interest.

As soon as we upload the paper on the arXiv (coming soon!), we will make available here the option to generate light cones and intensity maps for redshifts of interest, as well as analysis plots. Besides that, the user will be able to generate estimated luminosities on extended simulation boxes using Machine Learning techniques.

This algorithm was designed to be used on [SIMBA](http://simba.roe.ac.uk/) simulation datasets, but it will be expanded and generalized in the near future.

<br>

## developers and collaborators

*slick* was developed by **[Karolina Garcia](https://karolinagarcia.github.io/)** in collaboration with:
- Desika Narayanan
- Gergö Popping
- Anirudh Ravishankar
- Sagan Sutherland
- Melanie Kaasinen

<br>

## how it works

1) It takes a [SIMBA](http://simba.roe.ac.uk/) snapshot and its galaxy catalog (given by *[caesar](https://caesar.readthedocs.io/en/latest/)*), and writes the relevant information from all its particles in a pandas dataframe (which will be saved as "Basic_Characteristics*")

<br>

2) It creates a file with lists of cloud IDs (which will be saved as "Clouds_per_Core*"). Each line represents the IDs that will be passed to each core.

<br>

3) It runs radiative transfer equations on all the particles (or on a random set of those depending on the given parameters), and outputs a table with their CO, [CI] and [CII] luminosities.

<br>

## how to use
`setup.sh` marks the two frontend scripts as executable and adds the bin directory to the path.
This only needs to happen once ever.
```bash
bash setup.sh
source $HOME/.bashrc
```

To create a template slick project, run `slick_new.sh`. This creates [parameters.ini](#parameters-file) and `run.sh`.
`run.sh` calls slick_frontend.sh. It assumes slick_frontend.sh is on the path, so this will need to be changed if you haven't run `setup.sh` previously.
Running `sbatch run.sh` will enqueue the job for the slick_init step and the slick_run step (unless the skip_run option is set).

<br>

## options in the parameters file

The behavior of slick is configured by the parameters.ini file. The following describes the options currently available.
```ini
[snap]
; This is used when naming the clouds_per_core file.
boxsize=[int]
; The full path to the yt file being operated on.
ytfilename=[str]
; The full path to the caesar file being operated on.
caesarfilename=[str]

[sample]
; Either 'total' or 'randomized'.
; Determines whether slick should operate on all clouds or a sample of clouds.
mode=[str]
; How many galaxies should slick operate on.
; Only requireßd if mode = 'randomized'
n_galaxies_sample=[int]
; Minimum galaxy mass bound to operate on.
; Only required if mode = 'randomized'
min_mass=[float]
; Maximum galaxy mass bound to operate on.
; Only required if mode = 'randomized'
max_mass=[float]

[sbatch]
; This section takes any key value pairs which can be used in a job script as 
#SBATCH --key=value
; These are used to configure the slurm job for the slick_run step.
; The only parameter not configurable is array, which is set by slick to match the number of runs being prepared.
; The following is the default configuration
nodes=1
tasks-per-node=1
cpus-per-task=1
mem-per-cpu=8gb
time=96:00:00
output=/dev/null

[run]
; The directory which slick should output its files.
; This does not include any logs generated by slurm.
; To change the output of those, use the output parameter in [sbatch]
; Defaults to Output_Files
output_dir=[str]
; If true, only the slick_init step is run.
; The slick_run step can be triggered manually via `sbatch slick_run_jobscript.sh`
; Defaults to false
skip_run=[bool]
```
