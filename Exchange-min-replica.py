from imports import *
from input_file import *
#from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
#from definitions import *
from setup_pH_system import *
import subprocess
from pHrex import *


if restart != "ON":
    from createLambdaList import *


parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)

#if prep_replicas:
#    parallel_tempering.run_prep_replicas(MD_nsteps_replicas_prep, n_iter_replicas = 1000)
#if prep_lambdas:
#    parallel_tempering.run_prep_lambdas(MD_nsteps_lambdas_prep, n_iter_lambdas)

if HREMD == True:
    parallel_tempering.run_REX_EcpH(n_iter)
else:
    parallel_tempering.run()
