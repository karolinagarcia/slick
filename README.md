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

## Installing

SLICK depends on versions of DESPOTIC and CAESAR which are not available on pypi, so we'll have to install those ourselves.
Download these as well as SLICK itself.
```
  git clone git@github.com:karolinagarcia/slick.git
  git clone git@bitbucket.org:krumholz/despotic.git
  git clone git@github.com:dnarayanan/caesar.git
```

Make a python environment for slick to exist in:
```
  conda create -y --name slick python=3.10.4
  conda activate slick
```

DESPOTIC depends on libgsl 2.7. 
This can be resolved on HiPerGator by loading the corresponding module:
```
  module load gcc/12.2.0 gsl/2.7
```

We're going to install SLICK first.
It may seem weird to do this before the dependencies, but doing so in this order allows pip to install the dependencies that *are* on pypi (numpy, yt, etc.) for us.
```
  cd slick
  pip install .
  cd ..
```

Install DESPOTIC.
Doing so requires a patch to the makefile which allows the compilers to know where gsl is located on HiPerGator:
```
  cd despotic
  git checkout 182cd46d
  curl -L https://gist.githubusercontent.com/smsutherland/f12e6dac5bc91c5a227ea349dcce9098/raw/ | git apply
  python setup.py install
  cd ..
```

Install CAESAR:
```
  cd caesar
  git checkout da0dba1e
  python setup.py install
```

If all has gone well, you should be able to run ``slick -h`` and get a help message.

## how to use

Make sure you have the appropriate modules loaded and are in your slick conda environment:
```
  module load gcc/12.2.0 gsl/2.7
  conda activate slick
```

To initialize a slick run, use ``slick new``.
All the init step does is copy a ``parameters.ini`` and ``job.sh`` file into your current directory.
Presets are found in the ``slick/src/slick/presets/`` directory.
Currently only the default "narayanan" preset is shipped, but more can be made simply by adding "{name}.ini" and "{name}.sh" to the preset directory.
Any users not on HiPerGator are welcome to submit a PR adding their own presets.
``sbatch job.sh`` will queue the slick initialization step.
If the skip_run option is not set in the parameter file, the initialization step will automatically queue the run step of slick.

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

[module]
; This section takes any key value pairs which result in lmod modules being loaded at the beginning of the generated jobscript
gcc="12.2.0" ; Turns into module load gcc/12.2.0

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
