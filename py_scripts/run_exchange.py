from imports import *
from input_file import *
from setup_pH_system import *
import subprocess
from pHrex import *

parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)

parallel_tempering._mix_replicas()

exit
