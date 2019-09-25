from simtk import openmm, unit
from math import pi
import numpy as np

# Define VdW softcore parameters
switch_width=1*unit.angstroms
softcore_alpha = 0.5
softcore_a = 1
softcore_b = 1
softcore_c = 1
softcore_beta = 0.0
softcore_d = 1
softcore_e = 1
softcore_f = 1
# Define the replica exchange procedure
pH_low = 2
pH_high = 12
pH_step = 0.5
steps = 100
# Default pKa values
EpKa = 4.4
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
#Lambda lists
l_Glu = {}
l_His = {}
l_Lys = {}
l_Cys = {}
l_Asp = {}
l_Glu_special = {}
glu_char = {}
lys_char = {}
his_char_E = {}
his_char_D = {}
asp_char = {}
cys_char = {}
# GLU difference charges
glu_char["CG:"] = 0.07
glu_char["CD:"] = 0.13
glu_char["OE1:"] = 0.21
glu_char["OE2:"] = 0.15
# LYS difference charges
lys_char["CE:"] = 0.08
lys_char["HE1:"] = -0.025
lys_char["HE2:"] = -0.025
lys_char["NZ:"] = 0.66
lys_char["HZ1:"] = -0.01
lys_char["HZ2:"] = -0.01
# HSE difference charges
his_char_E["CB:"] = 0.03
his_char_E["ND1:"] = 0.19
#his_char_E["HD1:"] = 0.44
his_char_E["CG:"] = -0.03
his_char_E["CE1:"] = 0.07
his_char_E["HE1:"] = 0.05
his_char_E["NE2:"] = -0.15
his_char_E["HE2:"] = 0.12
his_char_E["CD2:"] = 0.24
his_char_E["HD2:"] = 0.04
# HSD difference charges
his_char_D["CB:"] = 0.04
his_char_D["ND1:"] = -0.15
his_char_D["HD1:"] = 0.12
his_char_D["CG:"] = 0.24
his_char_D["CE1:"] = 0.07
his_char_D["HE1:"] = 0.05
his_char_D["NE2:"] = 0.19
#his_char_D["HE2:"] = 0.44
his_char_D["CD2:"] = -0.03
his_char_D["HD2:"] = 0.03
# CYS difference charges
cys_char["CB:"] = 0.27
cys_char["SG:"] = 0.57
# ASP difference charges
asp_char["CB:"] = 0.07
asp_char["CG:"] = 0.13
asp_char["OD1:"] = 0.21
asp_char["OD2:"] = 0.15

for i in range(pH_low*10, pH_high*10, 5):
    p = 1-1/(1+10**(EpKa-i/10))
    l_Glu[i/10] = '%.4f'%p
    p = 1-1/(1+10**(HpKa-i/10))
    l_His[i/10] = '%.4f'%p
    p = 1-1/(1+10**(KpKa-i/10))
    l_Lys[i/10] = '%.4f'%p
    p = 1-1/(1+10**(CpKa-i/10))
    l_Cys[i/10] = '%.4f'%p
    p = 1-1/(1+10**(DpKa-i/10))
    l_Asp[i/10] = '%.4f'%p
pi = round(pi, 5)

pH_list = np.arange(pH_low, pH_high, 0.5)
