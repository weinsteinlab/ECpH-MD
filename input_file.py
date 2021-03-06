import os,sys
from imports import *
sys.path.append(os.path.join(sys.path[0],'py_scripts'))

# Job and Slurm settings
number_of_replicas=16     # Must be a multiple of # of GPUs per node
number_of_subjobs=4
jobName="example"         # No spaces
partitionName="edison"    # Slurm partition to run job on
number_of_GPUs_per_node=8 

# Input structure and coordinates #
psf_file = './inputFiles/bbl8.psf'
pdb_file = './inputFiles/bbl8.pdb'

# !!! Occurance of disulphide bonds #
disu = [66, 106, 119, 160]                       	# PDB residue numbers

# Output name #
output_name = 'test'
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

# pH range (range is inclusive for both low and high values)
pH_low = 8
pH_high = 9
pH_step = 1
replicas_per_pH = 8

# Special pKa values #
special_pKa_names = ['GLU89']
special_pKa_values = [7.3] 				# the same order and size as the names in special_pKa_names

# cpH setup #
minimize = True
n_min_steps = 1000

md_steps = 10000

