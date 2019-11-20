from imports import *
import subprocess, time
from definitions import *
from setup_pH_system import NBi, CNB, CB, pH_list
def manage_waters(pH_system_temp):
    nonbonded_force = pH_system_temp.getForces()[NBi]
    custom_electrostatics = [pH_system_temp.getForces()[CNB[2]], pH_system_temp.getForces()[CNB[3]]]
    net_charge = calc_system_charge(custom_electrostatics[0])
    if abs(net_charge) > 0.001:
        n_waters = int(abs(net_charge) / 0.001)
        print('Number of modified water oxygens ', n_waters)
        for w in waters:
            w = int(w)
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
            if charge._value != -0.834:
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
                nonbonded_force.setParticleParameters(w, wc * unit.elementary_charge, sigma, epsilon)
                for force in custom_electrostatics:
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(w)
                    force.setParticleParameters(w, [lambda_electrostatics, wc * unit.elementary_charge, sigma])

        if net_charge > 0:
            for w in r_waters:
                w = int(w)
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                wc = charge._value - 0.001
                nonbonded_force.setParticleParameters(w, wc * unit.elementary_charge, sigma, epsilon)
                for force in custom_electrostatics:
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(w)
                    force.setParticleParameters(w, [lambda_electrostatics, wc * unit.elementary_charge, sigma])


def pick_charges(res_list, residue, lambda_list):
    for resname in res_list:
        if str(resname) in residue:
            charge_list = res_list[resname]
        elif 'HSP' in residue:
            if float(lambda_list.at[(residue + '_sw', str(lambda_list.shape[1] - 1))]) == 1.0:
                charge_list = res_list['HSE']
            else:
                charge_list = res_list['HSD']

    return charge_list


def create_cpH_system(pH_system_temp, lambda_list):
    nonbonded_force = pH_system_temp.getForces()[NBi]
    for side_atom in side_atoms:
        side_atom = int(side_atom)
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
        residue = str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
        atom_name = str(psf.atom_list[side_atom].name)
        charge_list = pick_charges(res_list, residue, lambda_list)
        for atom in charge_list:
            if atom_name in atom:
                n_ch = charge._value - charge_list[atom] * (1 - float(lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                print(residue, '  ', atom_name, 'initial charge ', charge, ' assigned charge ', n_ch)

    print('CNB ', CNB)
    for f in CNB:
        force = pH_system_temp.getForces()[f]
        print('CNB force ', force)
        print('alchem protons in the loop ', alchem_protons)
        for proton in alchem_protons:
            proton = int(proton)
            residue = str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
            [lambda_sterics, sigma, epsilon] = force.getParticleParameters(proton)
            print('proton ', proton, ' lambda ', lambda_sterics, ' residue ', residue)
            if 'HSP' in residue:
                atom_name = str(psf.atom_list[proton].name)
                print('HSP atom name ', atom_name)
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
            print('proton index ', psf.atom_list[proton], ' assigned lambda value ', lambda_sterics)

        if force.getPerParticleParameterName(1) == 'charge':
            for side_atom in side_atoms:
                side_atom = int(side_atom)
                [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                residue = str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
                atom_name = str(psf.atom_list[side_atom].name)
                charge_list = pick_charges(res_list, residue, lambda_list)
                for atom in charge_list:
                    if atom_name in atom:
                        n_ch = charge - charge_list[atom] * (1 - float(lambda_list.at[(residue, str(lambda_list.shape[1] - 1))]))
                        force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                        print(residue, '  ', atom_name, 'initial charge ', charge, ' assigned charge ', n_ch)

    print('CB ', CB)
    for f in CB:
        force = pH_system_temp.getForces()[f]
        print('CB force ', force)
        for bond in range(force.getNumBonds()):
            [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
            if atom_i in alchem_protons:
                proton = atom_i
            elif atom_j in alchem_protons:
                proton = atom_j
            else:
                print('No alchemical protons in the bond') 
                break
            
            residue = str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
            atom_name = str(psf.atom_list[proton].name)
            if 'HSP' in residue:
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
            print('bond index ', bond, 'lambda value ', lambda_electrostatics)

        if force.getPerBondParameterName(1) == 'chargeprod':
            for bond in range(force.getNumBonds()):
                [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
                if atom_i in side_atoms:
                    charge_i = nonbonded_force.getParticleParameters(atom_i)[0]
                    charge_j = nonbonded_force.getParticleParameters(atom_j)[0]._value
                    force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, charge_i * charge_j, sigma])
                elif atom_j in side_atoms:
                    charge_i = nonbonded_force.getParticleParameters(atom_i)[0]._value
                    charge_j = nonbonded_force.getParticleParameters(atom_j)[0]
                    force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, charge_i * charge_j, sigma])


def create_cpH_exchange_system(pH_system, i, j, liex, lambda_list_i, ljex, lambda_list_j, new_name_i, new_name_j, proton_change):
    pH_system_temp = copy.deepcopy(pH_system)
    rep_ex = [i, j]
    for replica in rep_ex:
        nonbonded_force = pH_system_temp.getForces()[NBi]
        for side_atom in side_atoms:
            side_atom = int(side_atom)
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
            residue = str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
            if replica == i:
                if residue in proton_change:
                    print('Lambda swap will be attempted for ', str(residue), ' protonation state')
                    lambda_list = ljex
                    extended_lambda_list = lambda_list_j
                    new_name = new_name_j
                elif residue not in proton_change:
                    lambda_list = liex
                    extended_lambda_list = lambda_list_i
                    new_name = new_name_i
            elif replica == j:
                if residue in proton_change:
                    print('Lambda swap will be attempted for ', str(residue), ' protonation state')
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
            for proton in alchem_protons:
                proton = int(proton)
                residue = str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
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

                if 'HSP' in residue:
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
                for side_atom in side_atoms:
                    side_atom = int(side_atom)
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                    residue = str(psf.atom_list[side_atom].residue.resname) + str(psf.atom_list[side_atom].residue.idx)
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
                if atom_i in alchem_protons:
                    proton = atom_i
                elif atom_j in alchem_protons:
                    proton = atom_j
                residue = str(psf.atom_list[proton].residue.resname) + str(psf.atom_list[proton].residue.idx)
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
                if 'HSP' in residue:
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
                    if atom_i in side_atoms:
                        charge_i = nonbonded_force.getParticleParameters(atom_i)[0]
                        charge_j = nonbonded_force.getParticleParameters(atom_j)[0]._value
                        force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, charge_i * charge_j, sigma])
                    elif atom_j in side_atoms:
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

    def run_prep_lambdas(self, MD_nsteps_lambdas, n_iter_lambdas):
        self._MD_nsteps_lambdas = MD_nsteps_lambdas
        self._n_iter_lambdas = n_iter_lambdas
        for iteration in range(n_iter_lambdas):
            self._propagate_replicas(iteration, nsteps = MD_nsteps_lambdas)
            self._mix_lambdas()
    def run_prep_replicas(self, MD_nsteps_replicas, n_iter_replicas):
        self._MD_nsteps_replicas = MD_nsteps_replicas
        self._n_iter_replicas = n_iter_replicas
        for iteration in range(n_iter_replicas):
            self._propagate_replicas(iteration, nsteps = MD_nsteps_replicas)
            self._mix_replicas()
    def run(self, n_iter, MD_nsteps_lambdas, MD_nsteps_replicas):
        self._n_iter = n_iter
        self._MD_nsteps_lambdas = MD_nsteps_lambdas
        self._MD_nsteps_replicas = MD_nsteps_replicas
        for it in range(n_iter):
            self._propagate_replicas(it, MD_nsteps_replicas)
            print('\nbefore mix_replicas, after propagate\n')
            self._mix_replicas()
            self._propagate_replicas(it, MD_nsteps_lambdas)
            self._mix_lambdas()

    def _propagate_replicas(self, iteration, nsteps):
        slurmID = []
        for pH in self._pH_list:
            myCmd = 'sbatch submit_individual_replica.sh ' + str(pH) + ' ' + str(iteration) + ' ' + str(nsteps)
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
                    print('Sleepy time!\ni')
                    time.sleep(60)

    def _mix_lambdas(self, n_attempts=1):
        for attempt in range(n_attempts):
            i = randint(0, self._pH_list.size - 1)
            if i == 0:
                j = i + 1
            else:
                j = i - 1
            if n_residues_per_switch != None:
                number_atoms_change = round(len(list_alchem_residues) * n_residues_per_switch)
            else:
                number_atoms_change = round(len(list_alchem_residues) * 0.1)
          
            residues_change = random.sample(range(0, len(list_alchem_residues) - 1), number_atoms_change)
            proton_change = [None] * number_atoms_change
            for residue in range(len(residues_change)):
                proton_change[residue] = list_alchem_residues[residues_change[residue]]

            print('Exchange between replicas ', i, ' and ', j, ' protons ', proton_change)
            lambda_list_i = pd.read_csv(('lambda_list-' + str(self._pH_list[i]) + '.csv'), index_col=0)
            lambda_list_j = pd.read_csv(('lambda_list-' + str(self._pH_list[j]) + '.csv'), index_col=0)
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
            print('lambdas of replica i ', liex)
            print('lambdas of replica j ', ljex)
            pH_system_temp_i, pH_system_temp_j = create_cpH_exchange_system(self._pH_system, i, j, liex, lambda_list_i, ljex, lambda_list_j, new_name_i, new_name_j, proton_change)
            manage_waters(pH_system_temp)
            integrator_i = LangevinIntegrator(temperature, friction, dt)
            integrator_i.setConstraintTolerance(constraintTolerance)
            simulation_i = Simulation(topology, pH_system_temp_i, integrator_i, platform, platformProperties)
            simulation_i.loadState(str(output_name) + '-' + str(self._pH_list[i]) + '-state.xml')
            simulation_i.minimizeEnergy(maxIterations=1000)
            integrator_j = LangevinIntegrator(temperature, friction, dt)
            integrator_j.setConstraintTolerance(constraintTolerance)
            simulation_j = Simulation(topology, pH_system_temp_j, integrator_j, platform, platformProperties)
            simulation_j.loadState(str(output_name) + '-' + str(self._pH_list[j]) + '-state.xml')
            simulation_j.minimizeEnergy(maxIterations=1000)
            state_ji_lambda = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij_lambda = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
            energy_ji_lambda = state_ji_lambda.getKineticEnergy()._value + state_ji_lambda.getPotentialEnergy()._value
            energy_ij_lambda = state_ij_lambda.getKineticEnergy()._value + state_ij_lambda.getPotentialEnergy()._value
            energy_ii_lambda_file = pd.read_csv(str(output_name) + '-' + str(self._pH_list[i]) + '-energy-min.csv')
            energy_ii_lambda = float(energy_ii_lambda_file.iloc[:, -1])
            energy_jj_lambda_file = pd.read_csv(str(output_name) + '-' + str(self._pH_list[j]) + '-energy-min.csv')
            energy_jj_lambda = float(energy_jj_lambda_file.iloc[:, -1])
            print('Lambdas exchange energy ', energy_ii_lambda, energy_jj_lambda, energy_ij_lambda, energy_ji_lambda)
            log_p_accept_lambda = -(energy_ij_lambda + energy_ji_lambda) + energy_ii_lambda + energy_jj_lambda
            print(log_p_accept_lambda)
            if log_p_accept_lambda >= 0.0 or random.random() < math.exp(log_p_accept_lambda):
                print('Lambda exchange accepted')
                for residue in proton_change:
                    print('Shape ', lambda_list_i.shape[1], lambda_list_j.shape[1])
                    lambda_i = liex.at[(residue, str(lambda_list_i.shape[1]))]
                    lambda_j = ljex.at[(residue, str(lambda_list_j.shape[1]))]
                    liex.at[(residue, str(lambda_list_i.shape[1]))] = lambda_j
                    ljex.at[(residue, str(lambda_list_j.shape[1]))] = lambda_i

                df_output = pd.concat([lambda_list_i, liex], axis=1, sort=False)
                df_output.to_csv('lambda_list-' + str(self._pH_list[i]) + '.csv')
                df_output = pd.concat([lambda_list_j, ljex], axis=1, sort=False)
                df_output.to_csv('lambda_list-' + str(self._pH_list[j]) + '.csv')
            else:
                print('Lambda_exchange_rejected')

    def _mix_replicas(self, n_attempts=1):
        print("\nin mix_replicas")
        for attempt in range(n_attempts):
            i = randint(0, self._pH_list.size - 1)
            j = i
            while j == i:
                j = randint(0, self._pH_list.size - 1)

            rep_ex = [i, j]
            print('Exchange between replicas ', i, ' and ', j)
            lambda_list_i = pd.read_csv(('lambda_list-' + str(self._pH_list[i]) + '.csv'), index_col=0)
            lambda_list_j = pd.read_csv(('lambda_list-' + str(self._pH_list[j]) + '.csv'), index_col=0)
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
            print('lambdas of replica i ', liex)
            print('lambdas of replica j ', ljex)
            for replica in rep_ex:
                pH_system_temp = copy.deepcopy(self._pH_system)
                if replica == i:
                    lambda_list = lambda_list_i
                    new_name = new_name_i
                elif replica == j:
                    lambda_list = lambda_list_j
                    new_name = new_name_j
                else:
                    print('Error')
                create_cpH_system(pH_system_temp, lambda_list)
                manage_waters(pH_system_temp)
                if replica == i:
                    integrator_i = LangevinIntegrator(temperature, friction, dt)
                    integrator_i.setConstraintTolerance(constraintTolerance)
                    simulation_i = Simulation(topology, pH_system_temp, integrator_i, platform, platformProperties)
                    simulation_i.loadState(str(output_name) + '-' + str(self._pH_list[i]) + '-state.xml')
                    positions_i = simulation_i.context.getState(getPositions=True).getPositions()
                elif replica == j:
                    integrator_j = LangevinIntegrator(temperature, friction, dt)
                    integrator_j.setConstraintTolerance(constraintTolerance)
                    simulation_j = Simulation(topology, pH_system_temp, integrator_j, platform, platformProperties)
                    simulation_j.loadState(str(output_name) + '-' + str(self._pH_list[j]) + '-state.xml')
                    positions_j = simulation_j.context.getState(getPositions=True).getPositions()

            simulation_i.context.setPositions(positions_j)
            simulation_j.context.setPositions(positions_i)
            state_ji = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
            energy_ji = state_ji.getKineticEnergy()._value + state_ji.getPotentialEnergy()._value
            energy_ij = state_ij.getKineticEnergy()._value + state_ij.getPotentialEnergy()._value
            energy_ii_file = pd.read_csv(str(output_name) + '-' + str(self._pH_list[i]) + '-energy.csv')
            energy_ii = float(energy_ii_file.iloc[:, -1])
            energy_jj_file = pd.read_csv(str(output_name) + '-' + str(self._pH_list[j]) + '-energy.csv')
            energy_jj = float(energy_jj_file.iloc[:, -1])
            print('Replica exchange energy ', energy_ii, energy_jj, energy_ij, energy_ji)
            log_p_accept = -(energy_ij + energy_ji) + energy_ii + energy_jj
            if log_p_accept >= 0.0 or random.random() < math.exp(log_p_accept):
                print('Replica exchange accepted')
                simulation_i.context.setState(state_ij)
#                simulation_i.reporters.append(CheckpointReporter('test-' + str(self._pH_list[i]) + '.chk', 1))
                simulation_i.saveState(str(output_name) + '-' + str(self._pH_list[i]) + '-state.xml')
                simulation_i.context.createCheckpoint()
                simulation_j.context.setState(state_ji)
#                simulation_j.reporters.append(CheckpointReporter('test-' + str(self._pH_list[j]) + '.chk', 1))
                simulation_j.context.createCheckpoint()
                simulation_j.saveState(str(output_name) + '-' + str(self._pH_list[j]) + '-state.xml')
            else:
                print('Replica_exchange_rejected')
# okay decompiling pHrex_2.cpython-37.pyc