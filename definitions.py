from imports import *
from input_file import *
from setup_pH_system import *

params = CharmmParameterSet(top_file, par_file)

pH_list = np.arange(pH_low, pH_high, pH_step)

psf = CharmmPsfFile(psf_file)
psf.setBox(x_PBC_vector_length, y_PBC_vector_length, z_PBC_vector_length)

topology = psf.topology
EpKa = 4.4
EpKa2 = 7.3
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
l_Glu, l_His, l_Lys, l_Cys, l_Asp, l_Glu2, l_special = ({} for i in range(7))
glu_char, lys_char, his_char_E, his_char_D, asp_char, cys_char, res_list = ({} for i in range(7))
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

segment_list, alchem_residues, alchem_protons, list_alchem_residues, list_exchange_residues, side_atoms = ([] for i in range(6))
for remark_line in psf.title:
    if 'segment' in remark_line and 'NONE' not in remark_line:
        remark = remark_line.split(' ')
        for i in remark:
            if i =='segment':
                segment_list.append(remark[remark.index(i)+1])
print('Protein segments ', segment_list)

lys_atoms, cys_atoms, asp_atoms, glu_atoms, his_atoms_D, his_atoms_E = ([np.array([], dtype=(np.int16))]*len(segment_list) for i in range(6))
lys_side_atoms, cys_side_atoms, asp_side_atoms, glu_side_atoms, his_side_atoms = ([np.array([], dtype=(np.int16))]*len(segment_list) for i in range(5))

waters = np.array([], dtype=(np.int16))
CNB, CB = (np.array([], dtype=(np.int16)) for i in range(2))


with open(pdb_file) as template:
    for num, line in enumerate(template):
        column_line = line.split()
        for segment in segment_list:
            if segment in line:
                if 'LYS' in line:
                    if 'HZ3' in line:
                        lys_atoms[segment_list.index(segment)] = np.append(lys_atoms[segment_list.index(segment)], num - 1)
                    else:
                        for i in range(len(column_line)):
                            if lys_char.get(str(column_line[i])) != None:
                                lys_side_atoms[segment_list.index(segment)] = np.append(lys_side_atoms[segment_list.index(segment)], num - 1)
  
                elif 'GLU' in line:
                    if 'HE2' in line:
                        glu_atoms[segment_list.index(segment)] = np.append(glu_atoms[segment_list.index(segment)], num - 1)
                    else:
                        for i in range(len(column_line)):
                            if glu_char.get(str(column_line[i])) != None:
                                glu_side_atoms[segment_list.index(segment)] = np.append(glu_side_atoms[segment_list.index(segment)], num - 1)
  
                elif 'CYS' in line:
                    if 'HG1' in line:
                        cys_atoms[segment_list.index(segment)] = np.append(cys_atoms[segment_list.index(segment)], num - 1)
                    else:
                        for i in range(len(column_line)):
                            if cys_char.get(str(column_line[i])) != None:
                                if psf.atom_list[num - 1].residue.idx not in disu:
                                    cys_side_atoms[segment_list.index(segment)] = np.append(cys_side_atoms[segment_list.index(segment)], num - 1)
            
                elif 'ASP' in line:
                    if 'HD2' in line:
                        asp_atoms[segment_list.index(segment)] = np.append(asp_atoms[segment_list.index(segment)], num - 1)
                    else:
                        for i in range(len(column_line)):
                            if asp_char.get(str(column_line[i])) != None:
                                asp_side_atoms[segment_list.index(segment)] = np.append(asp_side_atoms[segment_list.index(segment)], num - 1)
  
                elif 'HSP' in line:
                    if 'HD1' in line:
                        his_atoms_E[segment_list.index(segment)] = np.append(his_atoms_E[segment_list.index(segment)], num - 1)
                    elif 'HE2' in line:
                        his_atoms_D[segment_list.index(segment)] = np.append(his_atoms_D[segment_list.index(segment)], num - 1)
                    else:
                        for i in range(len(column_line)):
                            if his_char_E.get(str(column_line[i])) != None or his_char_D.get(str(column_line[i])) != None:
                                his_side_atoms[segment_list.index(segment)] = np.append(his_side_atoms[segment_list.index(segment)], num - 1)
            elif 'TIP' in line:
                if 'OH2' in line:
                    waters = np.append(waters, num - 1)

for segment in range(len(segment_list)):

    alchem_residues.append(np.concatenate((lys_atoms[segment], his_atoms_E[segment], glu_atoms[segment], asp_atoms[segment], cys_atoms[segment])))
    alchem_protons.append(np.concatenate((lys_atoms[segment], his_atoms_E[segment], his_atoms_D[segment], glu_atoms[segment], asp_atoms[segment], cys_atoms[segment])))
    side_atoms.append(np.concatenate((lys_side_atoms[segment], his_side_atoms[segment], glu_side_atoms[segment], asp_side_atoms[segment], cys_side_atoms[segment])))
    alchem_residues[segment].sort()
    alchem_protons[segment].sort()
    side_atoms[segment].sort()

    for index in range(alchem_residues[segment].size):
        residue = segment_list[segment] + ': ' + str(psf.atom_list[(alchem_residues[segment][index] - 1)].residue.resname) + str(psf.atom_list[(alchem_residues[segment][index] - 1)].residue.idx)
        list_alchem_residues.append(residue)
        list_exchange_residues.append(residue)
        if 'HSP' in residue:
            list_exchange_residues.append(residue + str('_sw'))
print('Water oxygen ', waters.shape)


print('alchem_residues ', alchem_residues)

print('all alchemical proton indicies \n', alchem_protons)

print('side_atoms ', side_atoms)

print('List of residues with modified protonation state: \n ', list_exchange_residues)
print('Number of titratable sites: ', alchem_residues)

for i in np.arange(pH_low, pH_high, pH_step ):
    p = 1 - 1 / (1 + 10 ** (EpKa - i))
    l_Glu[i] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (HpKa - i))
    l_His[i] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (KpKa - i))
    l_Lys[i] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (CpKa - i))
    l_Cys[i] = '%.4f' % p
    p = 1 - 1 / (1 + 10 ** (DpKa - i))
    l_Asp[i] = '%.4f' % p

try:
    if len(special_pKa_names) == len(special_pKa_values):
        for n in range(len(special_pKa_names)):
            l_special[special_pKa_names[n]] = {}
            print(special_pKa_names[n])
            for pH in np.arange(pH_low, pH_high, pH_step):
                p = 1 - 1 / (1 + 10 ** (special_pKa_values[n] - pH))
                print(p)
                l_special[special_pKa_names[n]][pH] = '%.4f' % p
    else:
        print("Number of special residues and number of pKa values given do not match")
except NameError:
    print("Residues with special pKa values are not detected")

pi = round(pi, 5)

platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'0',  'Precision':'mixed'}
