from imports import *
from input_file import *
from setup_pH_system import *
#from setup_pH_system import pH_list

pH_list = np.arange(pH_low, pH_high, pH_step)

psf = CharmmPsfFile(psf_file)
psf.setBox(x_PBC_vector_length, y_PBC_vector_length, z_PBC_vector_length)

topology = psf.topology
print(topology.getUnitCellDimensions(), ' unit cell dimensions')

EpKa = 4.4
EpKa2 = 7.3
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
l_Glu = {}
l_His = {}
l_Lys = {}
l_Cys = {}
l_Asp = {}
l_Glu2 = {}
l_special = {}
glu_char = {}
lys_char = {}
his_char_E = {}
his_char_D = {}
asp_char = {}
cys_char = {}
res_list = {}
glu_char['CG'] = 0.07
glu_char['CD'] = 0.13
glu_char['OE1'] = 0.21
glu_char['OE2'] = 0.15
lys_char['CE'] = 0.08
lys_char['HE1'] = -0.025
lys_char['HE2'] = -0.025
lys_char['NZ'] = 0.66
lys_char['HZ1'] = -0.01
lys_char['HZ2'] = -0.01
his_char_E['CB'] = 0.03
his_char_E['ND1'] = 0.19
his_char_E['CG'] = -0.03
his_char_E['CE1'] = 0.07
his_char_E['HE1'] = 0.05
his_char_E['NE2'] = -0.15
his_char_E['HE2'] = 0.12
his_char_E['CD2'] = 0.24
his_char_E['HD2'] = 0.04
his_char_D['CB'] = 0.04
his_char_D['ND1'] = -0.15
his_char_D['HD1'] = 0.12
his_char_D['CG'] = 0.24
his_char_D['CE1'] = 0.07
his_char_D['HE1'] = 0.05
his_char_D['NE2'] = 0.19
his_char_D['CD2'] = -0.03
his_char_D['HD2'] = 0.03
cys_char['CB'] = 0.27
cys_char['SG'] = 0.57
asp_char['CB'] = 0.07
asp_char['CG'] = 0.13
asp_char['OD1'] = 0.21
asp_char['OD2'] = 0.15
res_list['LYS'] = lys_char
res_list['HSE'] = his_char_E
res_list['HSD'] = his_char_D
res_list['CYS'] = cys_char
res_list['ASP'] = asp_char
res_list['GLU'] = glu_char
lys_atoms = np.array([], dtype=(np.int16))
cys_atoms = np.array([], dtype=(np.int16))
asp_atoms = np.array([], dtype=(np.int16))
glu_atoms = np.array([], dtype=(np.int16))
his_atoms_D = np.array([], dtype=(np.int16))
his_atoms_E = np.array([], dtype=(np.int16))
lys_side_atoms = np.array([], dtype=(np.int16))
cys_side_atoms = np.array([], dtype=(np.int16))
glu_side_atoms = np.array([], dtype=(np.int16))
asp_side_atoms = np.array([], dtype=(np.int16))
his_side_atoms = np.array([], dtype=(np.int16))
waters = np.array([], dtype=(np.int16))
CNB = np.array([], dtype=(np.int16))
CB = np.array([], dtype=(np.int16))


for line in open(pdb_file):
    if 'LYS' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if 'HZ3' in column_line:
                if column_line[i] == 'HZ3':
                    lys_atoms = np.append(lys_atoms, int(column_line[(i - 2)]) - 1)
            elif lys_char.get(str(column_line[i])) != None:
                lys_side_atoms = np.append(lys_side_atoms, int(column_line[(i - 2)]) - 1)

    elif 'GLU' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if 'HE2' in column_line:
                if column_line[i] == 'HE2':
                    glu_atoms = np.append(glu_atoms, int(column_line[(i - 2)]) - 1)
            elif glu_char.get(str(column_line[i])) != None:
                glu_side_atoms = np.append(glu_side_atoms, int(column_line[(i - 2)]) - 1)

    elif 'CYS' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if 'HG1' in column_line:
                if column_line[i] == 'HG1':
                    cys_atoms = np.append(cys_atoms, int(column_line[(i - 2)]) - 1)
            if cys_char.get(str(column_line[i])) != None:
                if psf.atom_list[(int(column_line[(i - 2)]) - 1)].residue.idx not in disu:
                    cys_side_atoms = np.append(cys_side_atoms, int(column_line[(i - 2)]) - 1)

    elif 'ASP' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if 'HD2' in column_line:
                if column_line[i] == 'HD2':
                    asp_atoms = np.append(asp_atoms, int(column_line[(i - 2)]) - 1)
            elif asp_char.get(str(column_line[i])) != None:
                asp_side_atoms = np.append(asp_side_atoms, int(column_line[(i - 2)]) - 1)

    elif 'HSP' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if 'HD1' in column_line or 'HE2' in column_line:
                if column_line[i] == 'HD1':
                    his_atoms_E = np.append(his_atoms_E, int(column_line[(i - 2)]) - 1)
                elif column_line[i] == 'HE2':
                    his_atoms_D = np.append(his_atoms_D, int(column_line[(i - 2)]) - 1)
            elif his_char_E.get(str(column_line[i])) != None or his_char_D.get(str(column_line[i])) != None:
                his_side_atoms = np.append(his_side_atoms, int(column_line[(i - 2)]) - 1)

    elif 'TIP3' in line:
        column_line = line.split(' ')
        for i in range(len(column_line)):
            if column_line[i] == 'OH2':
                waters = np.append(waters, int(column_line[(i - 2)]) - 1)

#print('HIS ', his_atoms_E.shape)
print('HIS: ', his_atoms_D.shape)
print('GLU: ', glu_atoms.shape)
print('ASP: ', asp_atoms.shape)
print('CYS: ', cys_atoms.shape)
print('LYS: ', lys_atoms.shape)
print('LYS side: ', lys_side_atoms.shape)
print('CYS side: ', cys_side_atoms.shape)
print('GLU side: ', glu_side_atoms.shape)
print('ASP side: ', asp_side_atoms.shape)
print('HIS side: ', his_side_atoms.shape)
print('Water oxygen: ', waters.shape)
alchem_residues = np.concatenate((lys_atoms, his_atoms_E, glu_atoms, asp_atoms, cys_atoms))
alchem_residues = np.sort(alchem_residues)
alchem_protons = np.concatenate((lys_atoms, his_atoms_E, his_atoms_D, glu_atoms, asp_atoms, cys_atoms))
alchem_protons = np.sort(alchem_protons)
print('all alchemical proton indicies \n', alchem_protons)
list_alchem_residues = [None] * alchem_residues.size
side_atoms = np.concatenate((lys_side_atoms, his_side_atoms, glu_side_atoms, asp_side_atoms, cys_side_atoms))
side_atoms = np.sort(side_atoms)
for i in range(alchem_residues.size):
    list_alchem_residues[i] = str(psf.atom_list[(alchem_residues[i] - 1)].residue.resname) + str(psf.atom_list[(alchem_residues[i] - 1)].residue.idx)

print('List of residues with modified protonation state: \n', list_alchem_residues)
print('Number of titratable sites: ', alchem_residues.size)
for i in range(pH_low * 10, pH_high * 10, 5):
    p = 1 - 1 / (1 + 10 ** (EpKa - i / 10))
    l_Glu[i / 10] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (HpKa - i / 10))
    l_His[i / 10] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (KpKa - i / 10))
    l_Lys[i / 10] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (CpKa - i / 10))
    l_Cys[i / 10] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (DpKa - i / 10))
    l_Asp[i / 10] = '%.4f' % p

if len(special_pKa_names) == len(special_pKa_values) != 0:
    for n in range(len(special_pKa_names)):
        l_special[special_pKa_names[n]] = {}
        print(special_pKa_names[n])
        for pH in range(pH_low * 10, pH_high * 10, 5):
            p = 1 - 1 / (1 + 10 ** (special_pKa_values[n] - pH / 10))
            print(p)
            l_special[special_pKa_names[n]][pH / 10] = '%.4f' % p

pi = round(pi, 5)

