import sys

sys.path.insert(0, '../../Input_Files')
from create_cloudspercore_list import create_cloudspercore_list
from create_basic_characteristics_table import create_basic_table

from configparser import ConfigParser
config = ConfigParser()
config.read('parameters.ini')

date = config['date']['date']
boxsize = config['snap']['boxsize']
ytfilename = config['snap']['ytfilename']
caesarfilename = config['snap']['caesarfilename']
mode = config['sample']['mode']
n_galaxies_sample = config['sample']['n_galaxies_sample']
min_mass = config['sample']['min_mass']
max_mass = config['sample']['max_mass']

# Creates table with basic characteristics for all the clouds
create_basic_table(boxsize,ytfilename,caesarfilename,date)

###

# Creates table with basic characteristics for all the clouds (if mode='total'), or for a sample of them (mode = 'randomize')
create_cloudspercore_list(boxsize,ytfilename,caesarfilename,mode,n_galaxies_sample,min_mass,max_mass,date)