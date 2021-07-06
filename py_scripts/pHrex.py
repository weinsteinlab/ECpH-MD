import sys, os
sys.path.append('../')

CWD = os.getcwd()
sys.path.append(CWD)


from imports import *
import subprocess, time
#from input_file import pdb_state_files
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, calc_system_charge, create_force_particle, create_force_bond)
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
                print('NB forces: residue ', residue, ' proton name ', atom_name, ' lambda value ', lambda_sterics)

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


def create_cpH_exchange_system(pH_system, i, j, liex, lambda_list_i, ljex, lambda_list_j, new_name_i, new_name_j, proton_change):
    pH_system_temp = copy.deepcopy(pH_system)
    rep_ex = [i, j]
    for replica in rep_ex:
        nonbonded_force = pH_system_temp.getForces()[NBi]
        for segment in range(len(segment_list)):
            for side_atom in side_atoms[segment]:
                side_atom = int(side_atom)
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
                residue = segment_list[segment] + ': ' + str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
                if replica == i:
                    if residue in proton_change:
                        lambda_list = ljex
                        extended_lambda_list = lambda_list_j
                        new_name = new_name_j
                    elif residue not in proton_change:
                        lambda_list = liex
                        extended_lambda_list = lambda_list_i
                        new_name = new_name_i
                elif replica == j:
                    if residue in proton_change:
                        lambda_list = liex
                        extended_lambda_list = lambda_list_i
                        new_name = new_name_i
                    elif residue not in proton_change:
                        lambda_list = ljex
                        extended_lambda_list = lambda_list_j
                        new_name = new_name_j
 
                atom_name = str(psf.atom_list[side_atom].name)
                charge_list = pick_charges(res_list, residue, extended_lambda_list)
                for atom in charge_list:
                    if atom_name in atom:
                        n_ch = charge._value - charge_list[atom] * (1 - float(lambda_list.at[(residue, new_name)]))
                        nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)

        for f in CNB:
            force = pH_system_temp.getForces()[f]
            for segment in range(len(segment_list)):

                for proton in alchem_protons[segment]:
                    proton = int(proton)
                    residue = segment_list[segment] + ': ' + str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
                    if replica == i:
                        if residue in proton_change:
                            lambda_list = ljex
                            new_name = new_name_j
                        elif residue not in proton_change:
                            lambda_list = liex
                            new_name = new_name_i
                    elif residue == j:
                        if residue in proton_change:
                            lambda_list = liex
                            new_name = new_name_i
                        elif residue not in proton_change:
                            lambda_list = ljex
                            new_name = new_name_j
 
                    [lambda_value, parameter1, parameter2] = force.getParticleParameters(proton)
 
                    if 'HSP' in residue or 'HIS' in residue:
                        atom_name = str(psf.atom_list[proton].name)
                        if float(lambda_list.at[(residue + '_sw', new_name)]) == 0.0:
                            if 'HE2' in atom_name:
                                a = lambda_list.at[(residue, new_name)]
                            elif 'HD1' in atom_name:
                                a = 1.0
                        elif float(lambda_list.at[(residue + '_sw', new_name)]) == 1.0:
                            if 'HE2' in atom_name:
                                a = 1.0
                            elif 'HD1' in atom_name:
                                a = lambda_list.at[(residue, new_name)]
                    else:
                        lambda_value = lambda_list.at[(residue, new_name)]
 
                    force.setParticleParameters(proton, [lambda_value, parameter1, parameter2])

            if force.getPerParticleParameterName(1) == 'charge':
                for segment in range(len(segment_list)):
                    for side_atom in side_atoms[segment]:
                        side_atom = int(side_atom)
                        [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                        residue = segment_list[segment] + ': ' + str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
                        atom_name = str(psf.atom_list[side_atom].name)
                        if replica == i:
                            if residue in proton_change:
                                lambda_list = ljex
                                extended_lambda_list = lambda_list_j
                                new_name = new_name_j
                            elif residue not in proton_change:
                                lambda_list = liex
                                extended_lambda_list = lambda_list_i
                                new_name = new_name_i
                        elif replica == j:
                            if residue in proton_change:
                                lambda_list = liex
                                extended_lambda_list = lambda_list_i
                                new_name = new_name_i
                            elif residue not in proton_change:
                                lambda_list = ljex
                                extended_lambda_list = lambda_list_j
                                new_name = new_name_j
                                
                        charge_list = pick_charges(res_list, residue, extended_lambda_list)
                        for atom in charge_list:
                            if atom_name in atom:
                                n_ch = charge - charge_list[atom] * (1 - float(lambda_list.at[(residue, new_name)]))
                                force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])

        for f in CB:
            force = pH_system_temp.getForces()[f]
            for bond in range(force.getNumBonds()):
                [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
                for segment in range(len(segment_list)):
                    if atom_i in alchem_protons[segment]:
                        proton = atom_i
                        add = segment_list[segment] + ': '
                    elif atom_j in alchem_protons[segment]:
                        proton = atom_j
                        add = segment_list[segment] + ': '
                residue = add + str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
                atom_name = str(psf.atom_list[proton].name)
                if replica == i:
                    if residue in proton_change:
                        lambda_list = ljex
                        extended_lambda_list = lambda_list_j
                        new_name = new_name_j
                    elif residue not in proton_change:
                        lambda_list = liex
                        extended_lambda_list = lambda_list_i
                        new_name = new_name_i
                elif replica == j:
                    if residue in proton_change:
                        lambda_list = liex
                        extended_lambda_list = lambda_list_i
                        new_name = new_name_i
                    elif residue not in proton_change:
                        lambda_list = ljex
                        extended_lambda_list = lambda_list_j
                        new_name = new_name_j
                if 'HSP' in residue or 'HIS' in residue:
                    if float(lambda_list.at[(residue + '_sw', new_name)]) == 0.0:
                        if 'HE2' in atom_name:
                            lambda_electrostatics = lambda_list.at[(residue, new_name)]
                        elif 'HD1' in atom_name:
                            lambda_electrostatics = 1.0
                    elif float(lambda_list.at[(residue + '_sw', new_name)]) == 1.0:
                        if 'HE2' in atom_name:
                            lambda_electrostatics = 1.0
                        elif 'HD1' in atom_name:
                            lambda_electrostatics = lambda_list.at[(residue, new_name)]
                else:
                    lambda_electrostatics = lambda_list.at[(residue, new_name)]

                force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, chargeprod, sigma])

            if force.getPerBondParameterName(1) == 'chargeprod':
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

        if replica == i:
            pH_system_temp_i = pH_system_temp
        elif replica == j:
            pH_system_temp_j = pH_system_temp

    return (
     pH_system_temp_i, pH_system_temp_j)


class pHrex:

    def __init__(self, pH_system, pH_list):
        self._pH_system = pH_system
        self._pH_list = pH_list

    def run_REX_EcpH(self, n_iter):
        self._n_iter = n_iter
        for it in range(n_iter):
            self._propagate_replicas()
            self._mix_replicas()
    def run(self):
        self._propagate_replicas()

    def _propagate_replicas(self):
        slurmID = []
        for pH in self._pH_list:
            myCmd = 'sbatch -J ${SLURM_JOB_NAME}_${pH} submit_individual_replica.sh ' + str(pH)
            process = subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
            jobID = int(''.join(list(filter(str.isdigit, str(process.stdout)))))
            slurmID.append(jobID)

        for jobID in slurmID:
            while 1:
                myCmd = 'sacct -T -ojobid,state --noheader -j ' + str(jobID)
                process = subprocess.run(myCmd, shell=True, stdout=subprocess.PIPE)
                if 'FAILED' in str(process.stdout):
                    print('There was an error with Slurm JobID: ' + str(jobID))
                    quit()
                elif 'COMPLETED' in str(process.stdout):
                    print('Slurm JobID: ' + str(jobID) + ' COMPLETE')
                    break
                else:
                    print('Full process.stdout' + str(process.stdout) + '\n')
                    print('\nWaiting for Slurm JobID: '+str(jobID)+'\n')
                    time.sleep(60)


    def _mix_replicas(self, subjob_number, replicas_per_pH, number_of_replicas):
        
        for attempt in range(1):
            i = randint(0, self._pH_list.size - 1)
            j = i
            while j == i:
                j = randint(0, self._pH_list.size - 1)
          
            
            i_replica = randint(0, replicas_per_pH - 1)
            j_replica = randint(0, replicas_per_pH - 1)

            i_replica_directory = (i * replicas_per_pH) + i_replica
            j_replica_directory = (j * replicas_per_pH) + j_replica

            i_directory_path = './simulations/pH_' + str(self._pH_list[i]) + '_replica_number_' + str(i_replica_directory).zfill(4) + '/'
            j_directory_path = './simulations/pH_' + str(self._pH_list[j]) + '_replica_number_' + str(j_replica_directory).zfill(4) + '/'

            rep_ex = [i, j]
            print('\nSelected pHs for replica exchange: ' + str(self._pH_list[i]) + ' and ' + str(self._pH_list[j]) + '\n')
            print('\nSelected replica directories for replica exchange: ' + str(i_replica_directory) + ' and ' +  str(j_replica_directory) + '\n')
            lambda_list_i = pd.read_csv(('./lambdas/lambda_list-' + str(self._pH_list[i]) + '.csv'), index_col=0)
            lambda_list_j = pd.read_csv(('./lambdas/lambda_list-' + str(self._pH_list[j]) + '.csv'), index_col=0)
            lambda_iex = lambda_list_i.iloc[:, -1]
            lambda_jex = lambda_list_j.iloc[:, -1]
            lambda_iex = pd.DataFrame(lambda_iex)
            lambda_jex = pd.DataFrame(lambda_jex)
            name_i = str(lambda_iex.columns.tolist()[0])
            name_j = str(lambda_jex.columns.tolist()[0])
            new_name_i = str(lambda_list_i.shape[1])
            new_name_j = str(lambda_list_j.shape[1])
            liex = lambda_iex.rename(columns={name_i: new_name_i})
            ljex = lambda_jex.rename(columns={name_j: new_name_j})
            pH_system_temp = copy.deepcopy(self._pH_system)
            nonbonded_force = pH_system_temp.getForces()[7]
            for replica in rep_ex:
                pH_system_temp = copy.deepcopy(self._pH_system)
                if replica == i:
                    lambda_list = lambda_list_i
                    new_name = new_name_i
                elif replica == j:
                    lambda_list = lambda_list_j
                    new_name = new_name_j
                create_cpH_system(pH_system_temp, lambda_list)
                manage_waters(pH_system_temp)
                if replica == i:
                    integrator_i = LangevinIntegrator(temperature, friction, dt)
                    integrator_i.setConstraintTolerance(constraintTolerance)
                    platformProperties = {'DeviceIndex':'0',  'Precision':'mixed'}   
                    simulation_i = Simulation(topology, pH_system_temp, integrator_i, platform, platformProperties)
                    simulation_i.loadState(i_directory_path + str(output_name) + '-ph' + str(self._pH_list[i]) + '_replica_number_' + str(i_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-state.xml')
                    positions_i = simulation_i.context.getState(getPositions=True).getPositions()
                elif replica == j:
                    integrator_j = LangevinIntegrator(temperature, friction, dt)
                    integrator_j.setConstraintTolerance(constraintTolerance)
                    platformProperties = {'DeviceIndex':'1',  'Precision':'mixed'}  
                    simulation_j = Simulation(topology, pH_system_temp, integrator_j, platform, platformProperties)
                    #simulation_j.loadState(str(output_name) + '-' + str(self._pH_list[j]) + '-state.xml')
                    simulation_j.loadState(j_directory_path + str(output_name) + '-ph' + str(self._pH_list[j]) + '_replica_number_' + str(j_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-state.xml')
                    positions_j = simulation_j.context.getState(getPositions=True).getPositions()

            simulation_i.context.setPositions(positions_j)

            simulation_j.context.setPositions(positions_i)

            state_ji = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
    ### In current version of the code with OpenMM 7.5.0 the units of potential energy (kJ/mol) and Boltzmann constant (J/mol do not match
            energy_ji = state_ji.getPotentialEnergy()/(AVOGADRO_CONSTANT_NA*BOLTZMANN_CONSTANT_kB*temperature*1000)
            energy_ij = state_ij.getPotentialEnergy()/(AVOGADRO_CONSTANT_NA*BOLTZMANN_CONSTANT_kB*temperature*1000)
            energy_ii_file = pd.read_csv('./energies/' + str(output_name) + '-ph' + str(self._pH_list[i]) + '_replica_number_' + str(i_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-energy.csv')
            energy_ii = float(energy_ii_file.iloc[:, -1])
            #energy_jj_file = pd.read_csv(str(output_name) + '-' + str(self._pH_list[j]) + '-energy.csv')
            energy_jj_file = pd.read_csv('./energies/' + str(output_name) + '-ph' + str(self._pH_list[j]) + '_replica_number_' + str(j_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-energy.csv')
            energy_jj = float(energy_jj_file.iloc[:, -1])
            print('\nTotal energy of replica ', i_replica_directory, ': ', energy_ii, '\nTotal energy of replica ', j_replica_directory, ': ', energy_jj, '\nTotal energy of replica ', i_replica_directory, ' after replica exchange ', energy_ij, '\nTotal energy of replica ', j_replica_directory, ' after replica exchange ', energy_ji, '\n')
            log_p_accept = -(energy_ij + energy_ji) + energy_ii + energy_jj
            if log_p_accept >= 0.0 or random.random() < math.exp(log_p_accept):
                print('Replica exchange accepted')
                simulation_i.context.setState(state_ij)
                #simulation_i.saveState(str(output_name) + '-' + str(self._pH_list[i]) + '-state.xml')
                simulation_i.saveState(i_directory_path + str(output_name) + '-ph' + str(self._pH_list[i]) + '_replica_number_' + str(i_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-state.xml')
                simulation_i.context.createCheckpoint()
                simulation_j.context.setState(state_ji)
                simulation_j.context.createCheckpoint()
                #simulation_j.saveState(str(output_name) + '-' + str(self._pH_list[j]) + '-state.xml')
                simulation_j.saveState(j_directory_path + str(output_name) + '-ph' + str(self._pH_list[j]) + '_replica_number_' + str(j_replica_directory).zfill(4) + '-subjob' + str(subjob_number).zfill(4) + '-state.xml')
            else:
                print('Replica exchange rejected')
