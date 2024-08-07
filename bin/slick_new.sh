cat > parameters.ini << ENDOFFILE
[sample]
mode=total
#mode=randomize
#mode=single
#gal_ids=[8,16,28,74,78,89]

agn=False

#the following will only be applied when 'randomize' is on
#n_galaxies_sample=2
#min_mass=10**8
#max_mass=10**12
#min_mass=10**11
#max_mass=3*10**12

[snap]
boxsize=25
ytfilename=/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5
caesarfilename=/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5

[sbatch]
nodes=1
tasks-per-node=1
cpus-per-task=32
mem-per-cpu=8gb
time=96:00:00

[run]
n_zones = 16
max_cores = 92
# the calculation for the highest amount of cores is the total we wanna use divided by "cpus-per-task" in slick_run_jobscript file, example: 3000/32 = use 92
n_clouds_per_core = 120
basictable_dir=    #insert where basic tables should go to
output_dir=    #insert output directory here
skip_lumcalc=False
skip_basictable=False
#overwrite=False
#skip_ml=False
#if skip_ml = False, we need to use randomize to work
ENDOFFILE

cat > run.sh << ENDOFFILE
#!/bin/bash
#SBATCH --qos=narayanan-b
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=1:00:00

slick_frontend.sh parameters.ini
ENDOFFILE
