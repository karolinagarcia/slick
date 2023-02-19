cat > parameters.ini << ENDOFFILE
[sample]
mode=total

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
basictable_dir=    #insert where basic tables should go to
output_dir=    #insert output directory here
skip_lumcalc=False
skip_basictable=False
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
