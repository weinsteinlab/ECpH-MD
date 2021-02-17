from imports import *
from input_file import *
from setup_pH_system import *
import subprocess
from pHrex import *


if restart != "ON":
    from createLambdaList import *


parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)

if HREMD == True:
    parallel_tempering.run_REX_EcpH(n_iter)
else:
    parallel_tempering.run()
