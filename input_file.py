import sys
from imports import *
sys.path.append("/gpfs/u/home/SMMS/SMMSktsk/scratch/pH-Replica-Exchange/py_scripts/")

# Input structure and coordinates #
psf_file = './inputFiles/BB-pH6-wb-i.psf'
pdb_file = './inputFiles/BB-pH6-wb-i.pdb'

# !!! Occurance of disulphide bonds #
disu = [66, 106, 119, 160]                       	# PDB residue numbers

# Output name #
output_name = 'BBL-ecph'
dcdout_freq = 5000                               	# Number of steps between writing a dcd file

# Force field parameters #
top_file  ='./inputFiles/all_top.rtf'
par_file = './inputFiles/parameters.prm'

# Openmm MD simulation setup #
nonbondedMethod = PME
nonbondedCutoff = 12 * angstroms
switchDistance = 10 * angstroms
ewaldErrorTolerance = 0.0005

# Constraints #
constraints = HBonds
rigidWater = True
constraintTolerance = 1e-06

# Integration parameters #
dt = 0.002 * picoseconds
temperature = 310 * kelvin
friction = 1.0 / picosecond
integrator_type = "Langevin" 				# Verlet/Brownian/VariableVerlet/VariableLangevin

#err_tol = 0.001 # Error tolerance for Variable time step integrators
#thermostat_type = "Andersen" # Used with Verlet Integrator only

# Barostate type #
barostat_type = "MonteCarlo" 				# MonteCarloAnisotropic/MonteCarloMembrane
pressure = 1.01325*bar

# For MonteCarloAnisotropic Barostat 
# Define whether the dimensions fo periodic box are allowed to change size
#scale_X = False #
#scale_Y = False # True/False 
#scale_Z = False #

#barostat_freq = 100

#surface_tension = 0.0 					# the default surface tension acting on the system in bar*nm
#xymode = XYIsotropic 					# XYIsotropic / XYAnisotropic
#zmode = ZFree 						# ZFree / ZFixed / ConstantVolume

# Periodic cell size (nm) #
x_PBC_vector_length = 7.47
y_PBC_vector_length = 7.87
z_PBC_vector_length = 7.39

# pH range #
pH_low = 5
pH_high = 8 # exclusive--i.e. make this one larger than the last pH desired
pH_step = 1
replicas_per_pH = 2

# Special pKa values #
special_pKa_names = []
special_pKa_values = [] 				# the same order and size as the names in special_pKa_names

# Restart option #
restart = 'OFF' 					# ON/OFF


# cpH setup #
minimize = True
n_min_steps = 1000

md_steps = 40000

HREMD = False						# Steps within 1 MD cycle
