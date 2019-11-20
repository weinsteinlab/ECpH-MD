from imports import *
# Input structure and coordinates #
psf_file = './inputFiles/BB-pH6-wb-i.psf'
pdb_file = './inputFiles/BB-pH6-wb-i.pdb'
# !!! Occurance of disulphide bonds #
disu = [66, 106, 119, 160]
# Output name #
output_name = 'new_test'
dcdout_freq = 1000
# Force field parameters #
params = CharmmParameterSet('./inputFiles/all_top.rtf', './inputFiles/parameters.prm')
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
temperature = 298 * kelvin
friction = 1.0 / picosecond
# Periodic cell size #
x_PBC_vector_length = 7.47
y_PBC_vector_length = 7.87
z_PBC_vector_length = 7.39
# Cuda platform #
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'0',  'Precision':'mixed'}
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
pH_low = 2
pH_high = 11
pH_step = 0.5
# Default pKa values #
EpKa = 4.4
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
# Special pKa values #
special_pKa_names = ['GLU89']
special_pKa_values = [7.3]
# Restart option #
restart = 'OFF' # ON/OFF
# cpH setup #
n_min_steps = 5000
n_residues_per_switch = 0.1 # Number of residues for which lambda exchange attempt will be attempted
n_iter = 1000 # Number of iterations of MD - lambda-exchange - MD - replica-exchnage cycle
MD_nsteps_lambdas = 1000 # Number of MD steps before lambda-exchamge attempt
MD_nsteps_replicas = 10000 # Number of MD steps before replica-exchnage attempt

# Preparatory runs #
prep_replicas = False # True/False - run MD-REX cycle before the main cpH REX/lambda-EX cycle
#n_iter_replicas = 1000 # Number of iterations of the preliminary MD-REX cycle
#MD_nsteps_replicas_prep = 1000 # Number of MD steps for preliminary MD-REX cycle
prep_lambdas = False # True/Flase - run MD-lambda-EX cycle before the main cpH REX/lambda-EX cycle
#n_iter_lambdas = 1000 # Number of iterations of the preliminary MD-lambda-EX cycle
#MD_nsteps_lambdas_prep = 1000 # Number of MD steps for preliminary MD-l-EX cycle
