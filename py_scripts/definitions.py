from imports import *
from assign_prot_prob import * 
asp_atom_names = ['CB', 'CG', 'OD1', 'OD2']
cys_atom_names = ['CB', 'SG']
glu_atom_names = ['CG', 'CD', 'OE1', 'OE2']
hsd_atom_names = ['CB', 'ND1', 'HD1', 'CG', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2']
hse_atom_names = ['CB', 'ND1', 'CG', 'CE1', 'HE1', 'NE2', 'HE2', 'CD2', 'HD2']
lys_atom_names = ['CE', 'HE1', 'HE2', 'NZ', 'HZ1', 'HZ2']
atom_names = [asp_atom_names, cys_atom_names, glu_atom_names, hsd_atom_names, hse_atom_names, lys_atom_names]

asp_atom_delta_charges = [0.07, 0.13, 0.21, 0.15]
cys_atom_delta_charges = [0.27, 0.57]
glu_atom_delta_charges = [0.07, 0.13, 0.21, 0.15]
hsd_atom_delta_charges = [0.04, -0.15, 0.12, 0.24, 0.07, 0.05, 0.19, -0.03, 0.03]
hse_atom_delta_charges = [0.03, 0.19, -0.03, 0.07, 0.05, -0.15, 0.12, 0.24, 0.04]
lys_atom_delta_charges = [0.08, -0.025, -0.025, 0.66, -0.01, -0.01]
atom_delta_charges = [asp_atom_delta_charges, cys_atom_delta_charges, glu_atom_delta_charges, hsd_atom_delta_charges, hse_atom_delta_charges, lys_atom_delta_charges]

EpKa = 4.4
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
