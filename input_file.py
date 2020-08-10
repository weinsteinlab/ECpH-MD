from imports import *
# Input structure and coordinates #
psf_file = './inputFiles/BB-pH6-wb-i.psf'
pdb_file = './inputFiles/BB-pH6-wb-i.pdb'

# Load pH-states geometries #

# The number of elements in the list should be the same as the number of pH replicas/values.
# In this example 11 structures for pHs from 2.0 to 7.0 are in "closed" state (same pdb file)
# and 7 structures for pHs from 7.5 to 10.5 are in "opened" state (same pdb file)

read_state_geometries = False 				# True/False; Default False

pdb_state_files = ['./inputFiles/BB-pH6-wb-i.pdb']*11 +['./inputFiles/BB-pH6-8-wb-i.pdb']*7
                                                                                               
                                                                                               
# !!! Occurance of disulphide bonds #
disu = [66, 106, 119, 160]                       	# PDB residue numbers

# Output name #
output_name = 'bb-cph'
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

# Soft-core potential parameters for alchemical forces #
switch_width = 1 * unit.angstroms
softcore_alpha = 0.5
softcore_a = 1
softcore_b = 1
softcore_c = 1
softcore_beta = 0.0
softcore_d = 1
softcore_e = 1
softcore_f = 1

# pH range #
pH_low = 9
pH_high = 12
pH_step = 1

# Special pKa values #
special_pKa_names = ['GLU89']
special_pKa_values = [7.3] 				# the same order and size as the names in special_pKa_names

# Restart option #
restart = 'ON' 					# ON/OFF

# cpH setup #
n_min_steps = 1000
n_residues_per_switch = 0.4 				# Persentage of residues for which lambda exchange attempt will be attempted from 0 to 1 (default 0.1)
n_iter = 1000000 						# Number of iterations of MD - lambda-exchange - MD - replica-exchnage cycle

# Lambda exchnage parameters
MD_nsteps_lambdas = 50 				# Number of MD steps before lambda-exchamge attempt
n_attempts_lambdas = 1 					# Number of attempts for lambda-exchange

# Relica exchnage parameters
MD_nsteps_replicas = 50 				# Number of MD steps before replica-exchnage attempt
n_attempts_replicas = 1 				# Number of attempts for replica-exchange

# Preparatory runs #

# Prep replica-exchange
prep_replicas = False 					# True/False - run MD-REX cycle before the main cpH REX/lambda-EX cycle
n_iter_replicas = 1					# Number of iterations of the preliminary MD-REX cycle
MD_nsteps_replicas_prep = 60000000 			# Number of MD steps for preliminary MD-REX cycle

# Prep lambda-exchange
prep_lambdas = False 					# True/Flase - run MD-lambda-EX cycle before the main cpH REX/lambda-EX cycle
n_iter_lambdas = 1000 					# Number of iterations of the preliminary MD-lambda-EX cycle
MD_nsteps_lambdas_prep = 5000 				# Number of MD steps for preliminary MD-l-EX cycle
# ECpH run

MD_steps_cpH = 90000000
