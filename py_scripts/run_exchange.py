from imports import *
from input_file import *
from setup_pH_system import *
import subprocess
from pHrex import *

subjob_number = int(sys.argv[1])
replicas_per_pH = int(sys.argv[2])
number_of_replicas = int(sys.argv[3])

parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)

parallel_tempering._mix_replicas(subjob_number - 1, replicas_per_pH, number_of_replicas)

exit
