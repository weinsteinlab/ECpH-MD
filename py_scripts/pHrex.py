import sys, os
sys.path.append('../')

CWD = os.getcwd()
sys.path.append(CWD)


from imports import *
import subprocess, time
from fep_functions import calc_system_charge
from assign_prot_prob import *
from setup_pH_system import NBi, CNB, CB

def manage_waters(pH_system_temp):
    nonbonded_force = pH_system_temp.getForces()[NBi]
    custom_electrostatics = [pH_system_temp.getForces()[CNB[2]], pH_system_temp.getForces()[CNB[3]]]
    net_charge = calc_system_charge(custom_electrostatics[0])
    if abs(net_charge) > 0.001:
        print('\nNon-zero net charge. Modifying the charge of water oxygens.\n')
        n_waters = int(abs(net_charge) / 0.001)
        print('Number of modified water oxygens ', n_waters)
        for w in waters:
            w = int(w)
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
            if charge._value != -0.834:
                print('\nRestoring the initial charges for water oxygens\n')
                charge._value = -0.834
                nonbonded_force.setParticleParameters(w, charge, sigma, epsilon)
            for force in custom_electrostatics:
                [lambda_electrostatics, charge, sigma] = force.getParticleParameters(w)
                if charge != -0.834:
                    charge = -0.834
                    force.setParticleParameters(w, [lambda_electrostatics, charge, sigma])
        r_waters = list(np.random.choice(waters, n_waters))
        if net_charge < 0:
            for w in r_waters:
                w = int(w)
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                wc = charge._value + 0.001
                nonbonded_force.setParticleParameters(w, wc * elementary_charge, sigma, epsilon)
                for force in custom_electrostatics:
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(w)
                    force.setParticleParameters(w, [lambda_electrostatics, wc * elementary_charge, sigma])

        if net_charge > 0:
            for w in r_waters:
                w = int(w)
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                wc = charge._value - 0.001
                nonbonded_force.setParticleParameters(w, wc * elementary_charge, sigma, epsilon)
                for force in custom_electrostatics:
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(w)
                    force.setParticleParameters(w, [lambda_electrostatics, wc * elementary_charge, sigma])
        print('Final net charge : ', calc_system_charge(custom_electrostatics[0])) 


def pick_charges(res_list, residue, lambda_list):
    for resname in res_list:
        if str(resname) in residue:
            charge_list = res_list[resname]
        elif 'HSP' in residue or 'HIS' in residue:
            if float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 1.0:
                charge_list = res_list['HSE']
            else:
                charge_list = res_list['HSD']

    return charge_list


def create_cpH_system(pH_system_temp, lambda_list):
    nonbonded_force = pH_system_temp.getForces()[NBi]
    print('\nScaling the charges of atoms adjustent to alchemical protons for standard Nonbonded force\n')
    for segment in range(len(segment_list)):
    
        for side_atom in side_atoms[segment]:
            print(side_atom, ' ', segment)
            side_atom = int(side_atom)
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
            residue = segment_list[segment] + ': ' + str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
            atom_name = str(psf.atom_list[side_atom].name)
            charge_list = pick_charges(res_list, residue, lambda_list)
            for atom in charge_list:
                if atom_name in atom:
                    n_ch = charge._value - charge_list[atom] * (1 - float(lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]))
                    nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                    print(residue, '  ', atom_name, 'initial charge: ', charge, '; assigned scaled charge: ', n_ch)

    print('CNB ', CNB)
    for f in CNB:
        force = pH_system_temp.getForces()[f]

        print('\nScaling custom electrostatics and sterics nonbonded forces for alchemical protons\n') 
        for segment in range(len(segment_list)):
            for proton in alchem_protons[segment]:
                proton = int(proton)
                residue = segment_list[segment] + ': ' + str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
                [lambda_sterics, sigma, epsilon] = force.getParticleParameters(proton)
                atom_name = str(psf.atom_list[proton].name)
                if 'HSP' in residue or 'HIS' in residue:
                    if float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 0.0:
                        if 'HE2' in atom_name:
                            lambda_sterics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]
                        elif 'HD1' in atom_name:
                            lambda_sterics = 1.0
                    elif float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 1.0:
                        if 'HE2' in atom_name:
                            lambda_sterics = 1.0
                        elif 'HD1' in atom_name:
                            lambda_sterics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]
                else:
                    lambda_sterics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]
               
                force.setParticleParameters(proton, [lambda_sterics, sigma, epsilon])
                print('NB forces for alchemical protons after: residue ', residue, ' proton name ', atom_name, ' lambda value ', lambda_sterics)

        if force.getPerParticleParameterName(1) == 'charge':
            for segment in range(len(segment_list)):
                for side_atom in side_atoms[segment]:
                    side_atom = int(side_atom)
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                    residue = segment_list[segment] + ': ' + str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
                    atom_name = str(psf.atom_list[side_atom].name)
                    charge_list = pick_charges(res_list, residue, lambda_list)
                    for atom in charge_list:
                        if atom_name in atom:
                            n_ch = charge - charge_list[atom] * (1 - float(lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                            print(residue, '  ', atom_name, 'initial charge ', charge, ' assigned charge ', n_ch)

    for f in CB:
        force = pH_system_temp.getForces()[f]
        print('\nScaling custom electrostatic and steric bond forces for alchemical protons\n')
        for bond in range(force.getNumBonds()):
            [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
            for segment in range(len(segment_list)):
                if atom_i in alchem_protons[segment]:
                    proton = atom_i
                    add = segment_list[segment] + ': '
                elif atom_j in alchem_protons[segment]:
                    proton = atom_j
                    add = segment_list[segment] + ': '
                else:
                    print('No alchemical protons in the bond') 
                    continue
            
            residue = add + str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
            atom_name = str(psf.atom_list[proton].name)
            if 'HSP' in residue or 'HIS' in residue:
                if float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 0.0:
                    if 'HE2' in atom_name:
                        lambda_electrostatics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]
                    elif 'HD1' in atom_name:
                        lambda_electrostatics = 1.0
                elif float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 1.0:
                    if 'HE2' in atom_name:
                        lambda_electrostatics = 1.0
                    elif 'HD1' in atom_name:
                        lambda_electrostatics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]
            else:
                lambda_electrostatics = lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]

            force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, chargeprod, sigma])
            print('Bond forces: residue ', residue, ' proton name ', atom_name, ' lambda value ', lambda_electrostatics)

        if force.getPerBondParameterName(1) == 'chargeprod':
            print('\nScaling the charges of atoms adjustent to alchemical protons in Custom Bond electostatic force\n')
            for bond in range(force.getNumBonds()):
                [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
                for segment in range(len(segment_list)):
            
                    if atom_i in side_atoms[segment]:
                        charge_i = nonbonded_force.getParticleParameters(atom_i)[0]
                        charge_j = nonbonded_force.getParticleParameters(atom_j)[0]._value
                        force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, charge_i * charge_j, sigma])

                    elif atom_j in side_atoms[segment]:
                        charge_i = nonbonded_force.getParticleParameters(atom_i)[0]._value
                        charge_j = nonbonded_force.getParticleParameters(atom_j)[0]
                        force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, charge_i * charge_j, sigma])

                    else:
                        continue

