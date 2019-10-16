from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
from definitions import *
import subprocess

from createLambdaList import *

print(lambda_list)
print('Building system...')

from setup_pH_system import *
from pHrex_2 import *

parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)

#parallel_tempering.run(n_iter = 1000, nsteps = 50000)
parallel_tempering.run(n_iter = 1, nsteps = 200)

