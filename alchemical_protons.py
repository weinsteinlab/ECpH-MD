from imports import *

t = md.load('./inputFiles/fep-45-i01-1-wb-i.pdb')
top = t.topology

# Define alchemical protons
lys_atoms = top.select('resname LYS and name HZ3')
his_atoms_E = top.select('rescode H and name HD1')
his_atoms_D = top.select('rescode H and name HE2')
glu_atoms = top.select('resname GLU2 and name HE2')
asp_atoms = top.select('resname ASP2 and name HD2')
cys_atoms = top.select('resname CYS and name HG1')

# Select water molecules to keep the net charge
waters = t.top.select('water and type O')

#Define side atoms with flexible charge
lys_side_atoms = top.select('rescode K and name CE  HE2 HE3 NZ HZ1 HZ2')
his_side_atoms = top.select('rescode H and name CB ND1 CG CE1 HE1 NE2 CD2 HD2 HE2 HD1')
glu_side_atoms = top.select('resname GLU2 and name CG CD OE1 OE2')

side_atoms = np.concatenate((lys_side_atoms, his_side_atoms, glu_side_atoms))
side_atoms = np.sort(side_atoms)

all_atoms_HSE = np.concatenate((lys_atoms, his_atoms_E, glu_atoms))
all_atoms_HSE = np.sort(all_atoms_HSE)
all_atoms = np.concatenate((lys_atoms, his_atoms_E, his_atoms_D, glu_atoms))
all_atoms = np.sort(all_atoms)
list_alchem_residues = [None]*all_atoms_HSE.size

