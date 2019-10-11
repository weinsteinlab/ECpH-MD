from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
from definitions import *
import subprocess

from createLambdaList import *

print(lambda_list)
print('Building system...')

from setup_pH_system import *

class pHRex:

    def __init__(self, pH_system, pH_list):
        self._pH_system = pH_system
        self._pH_list = pH_list
        #self._n_iter = n_iter
#        self._nsteps = nsteps
        
    def run(self, n_iter, nsteps):
        self._n_iter = n_iter
        self._nsteps = nsteps
        for iteration in range(n_iter):
            self._propagate_replicas(iteration, nsteps)
            self._mix_lambdas()
            self._propagate_replicas(iteration, nsteps)
            self._mix_replicas()

    def _propagate_replicas(self, iteration, nsteps):
        for pH in self._pH_list:
            myCmd="python run_replica.py " + str(pH) + " " + str(iteration) + " " + str(nsteps)
            process = subprocess.Popen(myCmd, shell=True, stdout=subprocess.PIPE)
            process.wait()

    def _mix_lambdas(self, n_attempts=1):
        # Attempt to switch two replicas at random. 
        for attempt in range(n_attempts):
            # Select two replicas at random.
            i = randint(0, self._pH_list.size-1)
            if i == 0:
                j = i+1
            else:
                j = i-1
            rep_ex = [i, j]

            #Select 25% of alchemical proton for an exchange attempt

            number_atoms_change = round(len(list_alchem_residues)*0.25)
            residues_change = random.sample(range(0, len(list_alchem_residues)-1), number_atoms_change)
            proton_change = [None]*number_atoms_change
            for residue in range(len(residues_change)):
                proton_change[residue] = list_alchem_residues[residues_change[residue]]
            print('Exchange between replicas ', i, ' and ', j, ' protons ', proton_change)
            lambda_list_i = pd.read_csv('lambda_list-'+str(pH_list[i])+'.csv', index_col=0)
            lambda_list_j = pd.read_csv('lambda_list-'+str(pH_list[j])+'.csv', index_col=0)
            lambda_iex = lambda_list_i.iloc[:, -1] 
            lambda_jex = lambda_list_j.iloc[:, -1]
            lambda_iex = pd.DataFrame(lambda_iex)
            lambda_jex = pd.DataFrame(lambda_jex)
            name_i = str(lambda_iex.columns.tolist()[0])
            name_j = str(lambda_jex.columns.tolist()[0])
            new_name_i = str(lambda_list_i.shape[1])
            new_name_j = str(lambda_list_j.shape[1])
            liex = lambda_iex.rename(columns={name_i:new_name_i})
            ljex = lambda_jex.rename(columns={name_j:new_name_j})            
            pH_system_temp = copy.deepcopy(self._pH_system)            
            nonbonded_force = pH_system_temp.getForces()[8]
            print('lambdas of replica i ', liex)
            print('lambdas of replica j ', ljex)

            for replica in rep_ex:
                pH_system_temp = copy.deepcopy(self._pH_system)            
                nonbonded_force = pH_system_temp.getForces()[8]        
                for side_atom in side_atoms:
                    side_atom = int(side_atom)
                    [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
                    residue = str(top.atom(side_atom)).split('-')[0]
                    print('replica ', replica, ' residue ', residue)
                    if replica == i:
                        if residue in proton_change:
                            print ('Yes')
                            lambda_list = ljex
                            new_name = new_name_j
                        elif residue not in proton_change:
                            print ('No')
                            lambda_list = liex
                            new_name = new_name_i
                    elif replica == j:
                        if residue in proton_change:
                            print('Yes')
                            lambda_list = liex
                            new_name = new_name_i
                        elif residue not in proton_change:
                            print('No')
                            lambda_list = ljex
                            new_name = new_name_j
                    else:
                        print("Error")

                    atom_name = str(top.atom(side_atom)).split('-')[1]
                    if 'LYS' in residue:
                        for charge_delta in lys_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)               
                                print('LYS parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
                    if 'GLU' in residue:
                        for charge_delta in glu_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
                                print('GLU parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
                    if 'HIS' in residue:
                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                            his_char = his_char_E
                        else:
                            his_char = his_char_D
                        for charge_delta in his_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
                                print('HIS parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
                    if 'ASP' in residue:
                        for charge_delta in asp_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
                                print('ASP parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
                    if 'CYS' in residue:
                        for charge_delta in cys_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)       
                                print('CYS parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
                for f in range(9, 17):
                    force = pH_system_temp.getForces()[f]
                    if force.__class__.__name__ == 'CustomNonbondedForce':
                        # Assign labmdas to alchemical protons
                        for proton in (all_atoms):
                            proton = int(proton)
                            residue = str(top.atom(proton)).split('-')[0]
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
                            [a, b, c] = force.getParticleParameters(proton)
                            if 'HIS' in residue:
                                atom_name = str(top.atom(proton)).split('-')[1]
                                if float(lambda_list.at[residue+'_sw', new_name]) == 0.0:
                                    if 'HE2' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                                    elif 'HD1' in atom_name:
                                        a = 1.0
                                elif float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                    if 'HE2' in atom_name:
                                        a = 1.0
                                    elif 'HD1' in atom_name: 
                                        a = lambda_list.at[residue, new_name]
                            else:        
                                a = lambda_list.at[residue, new_name]

                            force.setParticleParameters(proton, [a, b, c])

                        # Assign transition charge to neighboring alchemical atoms
                        if force.getPerParticleParameterName(1) == 'charge':
                            for side_atom in (side_atoms):
                                side_atom = int(side_atom)
                                [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                                residue = str(top.atom(side_atom)).split('-')[0]
                                atom_name = str(top.atom(side_atom)).split('-')[1]
                                if replica == i:
                                    if residue in proton_change:
                                        lambda_list = ljex
                                        new_name = new_name_j
                                    elif residue not in proton_change:
                                        lambda_list = liex                            
                                        new_name = new_name_i
                                elif replica == j:
                                    if residue in proton_change:
                                        lambda_list = liex
                                        new_name = new_name_i
                                    elif residue not in proton_change:
                                        lambda_list = ljex
                                        new_name = new_name_j
                                if 'LYS' in residue:
                                    for charge_delta in lys_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                        
                                elif 'GLU' in residue:
                                    for charge_delta in glu_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                                elif 'HIS' in residue:
                                    if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                        his_char = his_char_E
                                    else:
                                        his_char = his_char_D
                                    for charge_delta in his_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                                elif 'ASP' in residue:
                                    for charge_delta in asp_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                                elif 'CYS' in residue:
                                    for charge_delta in cys_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                             
                    if force.__class__.__name__ == 'CustomBondForce':
                        # Assign labmdas to alchemical protons                                               
                        for bond in range(force.getNumBonds()):
                            [q, r, (a, b, c)] = force.getBondParameters(bond)
                            if q in all_atoms:
                                proton = q
                            elif r in all_atoms:
                                proton = r
                            residue = str(top.atom(proton)).split('-')[0]
                            atom_name = str(top.atom(proton)).split('-')[1]
                            if replica == i:
                                if residue in proton_change:
                                    lambda_list = ljex
                                    new_name = new_name_j
                                elif residue not in proton_change:
                                    lambda_list = liex
                                    new_name = new_name_i
                            elif replica == j:
                                if residue in proton_change:
                                    lambda_list = liex
                                    new_name = new_name_i
                                elif residue not in proton_change:
                                    lambda_list = ljex
                                    new_name = new_name_j
                            if 'HIS' in residue:
                                if float(lambda_list.at[residue+'_sw', new_name]) == 0.0:
                                    if 'HE2' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                                    elif 'HD1' in atom_name:
                                        a = 1.0
                                elif float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                    if 'HE2' in atom_name:
                                        a = 1.0
                                    elif 'HD1' in atom_name:
                                        a = lambda_list.at[residue, new_name]                          
                            else:
                                a = lambda_list.at[residue, new_name]

                            force.setBondParameters(bond, q, r, [a, b, c])                                               

                        # Assign transition charge to neighboring alchemical atoms
                        if force.getPerBondParameterName(1) == 'chargeprod':
                            for bond in range(force.getNumBonds()):
                                [q, r, (a, b, c)] = force.getBondParameters(bond)
                                #b = round(b, 5)
                                if q in side_atoms:
                                    side_atom = q
                                    residue = str(top.atom(side_atom)).split('-')[0]
                                    atom_name = str(top.atom(side_atom)).split('-')[1]
                                    if replica == i: 
                                        if residue in proton_change:
                                            lambda_list = ljex
                                            new_name = new_name_j
                                        elif residue not in proton_change:
                                            lambda_list = liex                              
                                            new_name = new_name_i
                                    elif replica == j:
                                        if residue in proton_change:
                                            lambda_list = liex
                                            new_name = new_name_i
                                        elif residue not in proton_change:
                                            lambda_list = ljex
                                            new_name = new_name_j
                                    if 'LYS' in residue:
                                        for charge_delta in lys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])                                              
                                    if 'GLU' in residue:
                                        for charge_delta in glu_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'HIS' in residue:
                                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                            his_char = his_char_E
                                        else:
                                            his_char = his_char_D
                                        for charge_delta in his_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'ASP' in residue:
                                        for charge_delta in asp_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'CYS' in residue:
                                        for charge_delta in cys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c]) 
                                elif r in side_atoms:
                                    side_atom = r
                                    residue = str(top.atom(side_atom)).split('-')[0]
                                    atom_name = str(top.atom(side_atom)).split('-')[1]
                                    if replica == i:
                                        if residue in proton_change:
                                            lambda_list = ljex
                                            new_name = new_name_j
                                        elif residue not in proton_change:
                                            lambda_list = liex
                                            new_name = new_name_i
                                    elif replica == j:
                                        if residue in proton_change:
                                            lambda_list = liex
                                            new_name = new_name_i
                                        elif residue not in proton_change:
                                            lambda_list = ljex
                                            new_name = new_name_j

                                    if 'LYS' in residue:
                                        for charge_delta in lys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'GLU' in residue:
                                        for charge_delta in glu_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'HIS' in residue:
                                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                            his_char = his_char_E
                                        else:
                                            his_char = his_char_D
                                        for charge_delta in his_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'ASP' in residue:
                                        for charge_delta in asp_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'CYS' in residue:
                                        for charge_delta in cys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])

                el_nb = [pH_system_temp.getForces()[11], pH_system_temp.getForces()[12]]
                net_charge = calc_system_charge(el_nb[0])
                if abs(net_charge) > 0.001:
                    n_waters = int(abs(net_charge)/0.001)
                    print('Number of modified water oxygens ', n_waters)
                    for w in waters:
                        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(int(w))
                        if charge._value != -0.834:
                            charge._value = -0.834
                            nonbonded_force.setParticleParameters(i, charge, sigma, epsilon)
                        for force in el_nb:
                            [l, charge, sigma] = force.getParticleParameters(int(w))
                            if charge != -0.834:
                                charge = -0.834
                                force.setParticleParameters(i, [l, charge, sigma])
                    r_waters = list(np.random.choice(waters, n_waters))
                    if net_charge < 0:
                        for w in r_waters:
                            w = int(w)            
                            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                            wc = charge._value + 0.001
                            nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
                            for force in el_nb:
                                [l, charge, sigma] = force.getParticleParameters(w)
                                force.setParticleParameters(w, [l, wc*unit.elementary_charge, sigma])
                    if net_charge > 0:
                        for w in r_waters:
                            w = int(w)            
                            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                            wc = charge._value - 0.001        
                            nonbonded_force.setParticleParameters(w, wc*unit.elementary_charge, sigma, epsilon)
                            for force in el_nb:
                                [l, charge, sigma] = force.getParticleParameters(w)
                                force.setParticleParameters(w, l, wc*unit.elementary_charge, sigma)                    

            #Calculate total energy for selected replicas
            
                if replica == i:
                    
                    integrator_i = LangevinIntegrator(temperature, friction, dt)
                    integrator_i.setConstraintTolerance(constraintTolerance)
                    simulation_i = Simulation(topology, pH_system_temp, integrator_i, platform, platformProperties)
                    simulation_i.loadState('test-'+str(self._pH_list[i])+'.xml')
                    simulation_i.minimizeEnergy(maxIterations = 1000)
#                    positions_i = simulation_i.context.getState(getPositions=True).getPositions()

                elif replica == j:
                    
                    integrator_j = LangevinIntegrator(temperature, friction, dt)
                    integrator_j.setConstraintTolerance(constraintTolerance)
                    simulation_j = Simulation(topology, pH_system_temp, integrator_j, platform, platformProperties)
                    simulation_j.loadState('test-'+str(self._pH_list[j])+'.xml')
                    simulation_j.minimizeEnergy(maxIterations=1000)
#                    positions_j = simulation_j.context.getState(getPositions=True).getPositions()
#            simulatiion_i.context.setPositions(positions_j)
#            simulation_j.context.setPositions(positions_i)
            state_ji_lambda = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij_lambda = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            energy_ji_lambda = state_ji_lambda.getKineticEnergy()._value + state_ji_lambda.getPotentialEnergy()._value
            energy_ij_lambda = state_ij_lambda.getKineticEnergy()._value + state_ij_lambda.getPotentialEnergy()._value
            #Calculate total energy for replicas with exchnaged lambda values
            energy_ii_lambda_file = pd.read_csv('energy-min-'+str(self._pH_list[i])+'.csv')
            energy_ii_lambda = float(energy_ii_lambda_file.iloc[:, -1])
            energy_jj_lambda_file = pd.read_csv('energy-min-'+str(self._pH_list[j])+'.csv')
            energy_jj_lambda = float(energy_jj_lambda_file.iloc[:, -1])
            print("Lambdas exchange energy ", energy_ii_lambda, energy_jj_lambda, energy_ij_lambda, energy_ji_lambda) 
            # Accept or reject the swap.
            log_p_accept_lambda = - (energy_ij_lambda + energy_ji_lambda) + energy_ii_lambda + energy_jj_lambda
            print(log_p_accept_lambda)
            if log_p_accept_lambda >= 0.0 or random.random() < math.exp(log_p_accept_lambda):
                print('Lambda exchange accepted')
                # Swap states in replica slots i and j.
#                simulation_i.context.setState(state_ij)
#                simulation_i.reporters.append(CheckpointReporter('test-'+str(self._pH_list[i])+'.chk', 1))
#                simulation_i.saveState('test-'+str(self._pH_list[i])+'.xml')
#                simulation_i.context.createCheckpoint()
#                simulation_j.context.setState(state_ji)
#                simulation_j.reporters.append(CheckpointReporter('test-'+str(self._pH_list[j])+'.chk', 1))
#                simulation_j.context.createCheckpoint()
#                simulation_j.saveState('test-'+str(self._pH_list[j])+'.xml')
                for residue in proton_change:
#                    proton = int(proton)
#                    residue = str(top.atom(proton)).split('-')[0]
                    print('Shape ', lambda_list_i.shape[1], lambda_list_j.shape[1])
                    lambda_i = liex.at[residue, str(lambda_list_i.shape[1])]
                    lambda_j = ljex.at[residue, str(lambda_list_j.shape[1])]
                    liex.at[residue, str(lambda_list_i.shape[1])] = lambda_j
                    ljex.at[residue, str(lambda_list_j.shape[1])] = lambda_i
#                    df = pd.DataFrame.from_dict(lambda_list_i, orient='index', columns=[new_name])
#                    df_input = pd.read_csv('lambda_list-'+str(self._pH_list[i])+'.csv')
                df_output = pd.concat([lambda_list_i, liex], axis = 1, sort = False)
                df_output.to_csv('lambda_list-'+str(self._pH_list[i])+'.csv')
#                    df = pd.DataFrame.from_dict(lambda_list_j, orient='index', columns=[new_name])
#                    df_input = pd.read_csv('lambda_list-'+str(self._pH_list[j])+'.csv')
                df_output = pd.concat([lambda_list_j, ljex], axis = 1, sort = False)
                df_output.to_csv('lambda_list-'+str(self._pH_list[j])+'.csv')
            else:
                print('Lambda_exchange_rejected')


    def _mix_replicas(self, n_attempts=1):
        # Attempt to switch two replicas at random. 
        for attempt in range(n_attempts):
            # Select two replicas at random.
            i = randint(0, self._pH_list.size-1)
            j = i
            while j == i:
                j = randint(0, self._pH_list.size-1)
            rep_ex = [i, j]
            #Select 25% of alchemical proton for an exchange attempt
#            number_atoms_change = round(len(list_alchem_residues)*0.25)
#            residues_change = random.sample(range(0, len(list_alchem_residues)-1), number_atoms_change)
#            proton_change = [None]*number_atoms_change
#            for residue in range(len(residues_change)):
#                proton_change[residue] = list_alchem_residues[residues_change[residue]]
            print('Exchange between replicas ', i, ' and ', j)
            lambda_list_i = pd.read_csv('lambda_list-'+str(pH_list[i])+'.csv', index_col=0)
            lambda_list_j = pd.read_csv('lambda_list-'+str(pH_list[j])+'.csv', index_col=0)
            lambda_iex = lambda_list_i.iloc[:, -1]
            lambda_jex = lambda_list_j.iloc[:, -1]
            lambda_iex = pd.DataFrame(lambda_iex)
            lambda_jex = pd.DataFrame(lambda_jex)
            name_i = str(lambda_iex.columns.tolist()[0])
            name_j = str(lambda_jex.columns.tolist()[0])
            new_name_i = str(lambda_list_i.shape[1])
            new_name_j = str(lambda_list_j.shape[1])
            liex = lambda_iex.rename(columns={name_i:new_name_i})
            ljex = lambda_jex.rename(columns={name_j:new_name_j})
            pH_system_temp = copy.deepcopy(self._pH_system)
            nonbonded_force = pH_system_temp.getForces()[8]
            print('lambdas of replica i ', liex)
            print('lambdas of replica j ', ljex)

            for replica in rep_ex:
                pH_system_temp = copy.deepcopy(self._pH_system)
                nonbonded_force = pH_system_temp.getForces()[8]
                for side_atom in side_atoms:
                    side_atom = int(side_atom)
                    [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
                    residue = str(top.atom(side_atom)).split('-')[0]
                    print('replica ', replica, ' residue ', residue)
                    if replica == i:
                        lambda_list = liex
                        new_name = new_name_i
                    elif replica == j:
                        lambda_list = ljex
                        new_name = new_name_j
                    else:
                        print("Error")

                    atom_name = str(top.atom(side_atom)).split('-')[1]
                    if 'LYS' in residue:
                        for charge_delta in lys_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                                print('LYS parameters ', side_atom, n_ch, ' old charge: ', charge._value)
                    if 'GLU' in residue:
                        for charge_delta in glu_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                                print('GLU parameters ', side_atom, n_ch, ' old charge: ', charge._value)
                    if 'HIS' in residue:
                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                            his_char = his_char_E
                        else:
                            his_char = his_char_D
                        for charge_delta in his_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                                print('HIS parameters ', side_atom, n_ch, ' old charge: ', charge._value)
                    if 'ASP' in residue:
                        for charge_delta in asp_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                                print('ASP parameters ', side_atom, n_ch, ' old charge: ', charge._value)
                    if 'CYS' in residue:
                        for charge_delta in cys_char:
                            if atom_name in charge_delta:
                                n_ch = charge._value - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
                                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
                                print('CYS parameters ', side_atom, n_ch, ' old charge: ', charge._value)
                for f in range(9, 17):
                    force = pH_system_temp.getForces()[f]
                    if force.__class__.__name__ == 'CustomNonbondedForce':
                        # Assign labmdas to alchemical protons
                        for proton in (all_atoms):
                            proton = int(proton)
                            residue = str(top.atom(proton)).split('-')[0]
                            if replica == i:
                                lambda_list = liex
                                new_name = new_name_i
                            elif residue == j:
                                lambda_list = ljex
                                new_name = new_name_j
                            [a, b, c] = force.getParticleParameters(proton)
                            if 'HIS' in residue:
                                atom_name = str(top.atom(proton)).split('-')[1]
                                if float(lambda_list.at[residue+'_sw', new_name]) == 0.0:
                                    if 'HE2' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                                    elif 'HD1' in atom_name:
                                        a = 1.0
                                elif float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                    if 'HE2' in atom_name:
                                        a = 1.0
                                    elif 'HD1' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                            else:
                                a = lambda_list.at[residue, new_name]

                            force.setParticleParameters(proton, [a, b, c])

                        # Assign transition charge to neighboring alchemical atoms
                        if force.getPerParticleParameterName(1) == 'charge':
                            for side_atom in (side_atoms):
                                side_atom = int(side_atom)
                                [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                                residue = str(top.atom(side_atom)).split('-')[0]
                                atom_name = str(top.atom(side_atom)).split('-')[1]
                                if replica == i:
                                    lambda_list = liex
                                    new_name = new_name_i
                                elif replica == j:
                                    lambda_list = ljex
                                    new_name = new_name_j
                                if 'LYS' in residue:
                                    for charge_delta in lys_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                        print(charge, ' ', n_ch)
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                                elif 'GLU' in residue:
                                    for charge_delta in glu_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                        print(charge, ' ', n_ch)
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                                elif 'HIS' in residue:
                                    if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                        his_char = his_char_E
                                    else:
                                        his_char = his_char_D
                                    for charge_delta in his_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                        print(charge, ' ', n_ch)
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                                elif 'ASP' in residue:
                                    for charge_delta in asp_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                        print(charge, ' ', n_ch)
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                                elif 'CYS' in residue:
                                    for charge_delta in cys_char:
                                        if atom_name in charge_delta:
                                            n_ch = charge - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                        print(charge, ' ', n_ch)
                                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])
                    if force.__class__.__name__ == 'CustomBondForce':
                        # Assign labmdas to alchemical protons                                               
                        for bond in range(force.getNumBonds()):
                            [q, r, (a, b, c)] = force.getBondParameters(bond)
                            if q in all_atoms:
                                proton = q
                            elif r in all_atoms:
                                proton = r
                            residue = str(top.atom(proton)).split('-')[0]
                            atom_name = str(top.atom(proton)).split('-')[1]
                            if replica == i:
                                lambda_list = liex
                                new_name = new_name_i
                            elif replica == j:
                                lambda_list = ljex
                                new_name = new_name_j
                            if 'HIS' in residue:
                                if float(lambda_list.at[residue+'_sw', new_name]) == 0.0:
                                    if 'HE2' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                                    elif 'HD1' in atom_name:
                                        a = 1.0
                                elif float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                    if 'HE2' in atom_name:
                                        a = 1.0
                                    elif 'HD1' in atom_name:
                                        a = lambda_list.at[residue, new_name]
                            else:
                                a = lambda_list.at[residue, new_name]

                            force.setBondParameters(bond, q, r, [a, b, c])

                        # Assign transition charge to neighboring alchemical atoms
                        if force.getPerBondParameterName(1) == 'chargeprod':
                            for bond in range(force.getNumBonds()):
                                [q, r, (a, b, c)] = force.getBondParameters(bond)
                                #b = round(b, 5)
                                if q in side_atoms:
                                    side_atom = q
                                    residue = str(top.atom(side_atom)).split('-')[0]
                                    atom_name = str(top.atom(side_atom)).split('-')[1]
                                    if replica == i:
                                        lambda_list = liex
                                        new_name = new_name_i
                                    elif replica == j:
                                        lambda_list = ljex
                                        new_name = new_name_j
                                    if 'LYS' in residue:
                                        for charge_delta in lys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'GLU' in residue:
                                        for charge_delta in glu_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'HIS' in residue:
                                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                            his_char = his_char_E
                                        else:
                                            his_char = his_char_D
                                        for charge_delta in his_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'ASP' in residue:
                                        for charge_delta in asp_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'CYS' in residue:
                                        for charge_delta in cys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                elif r in side_atoms:
                                    side_atom = r
                                    residue = str(top.atom(side_atom)).split('-')[0]
                                    atom_name = str(top.atom(side_atom)).split('-')[1]
                                    if replica == i:
                                        lambda_list = liex
                                        new_name = new_name_i
                                    elif replica == j:
                                        lambda_list = ljex
                                        new_name = new_name_j

                                    if 'LYS' in residue:
                                        for charge_delta in lys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'GLU' in residue:
                                        for charge_delta in glu_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'HIS' in residue:
                                        if float(lambda_list.at[residue+'_sw', new_name]) == 1.0:
                                            his_char = his_char_E
                                        else:
                                            his_char = his_char_D
                                        for charge_delta in his_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'ASP' in residue:
                                        for charge_delta in asp_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                                    if 'CYS' in residue:
                                        for charge_delta in cys_char:
                                            if atom_name in charge_delta:
                                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, new_name]))
    #                                            print(b, ' ', n_ch)
                                                force.setBondParameters(bond, q, r, [a, n_ch, c])

                el_nb = [pH_system_temp.getForces()[11], pH_system_temp.getForces()[12]]
                net_charge = calc_system_charge(el_nb[0])
                if abs(net_charge) > 0.001:
                    n_waters = int(abs(net_charge)/0.001)
                    print('Number of modified water oxygens ', n_waters)
                    for w in waters:
                        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(int(w))
                        if charge._value != -0.834:
                            charge._value = -0.834
                            nonbonded_force.setParticleParameters(i, charge, sigma, epsilon)
                        for force in el_nb:
                            [l, charge, sigma] = force.getParticleParameters(int(w))
                            if charge != -0.834:
                                charge = -0.834
                                force.setParticleParameters(i, [l, charge, sigma])
                    r_waters = list(np.random.choice(waters, n_waters))
                    if net_charge < 0:
                        for w in r_waters:
                            w = int(w)
                            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                            wc = charge._value + 0.001
                            nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
                            for force in el_nb:
                                [l, charge, sigma] = force.getParticleParameters(w)
                                force.setParticleParameters(w, [l, wc*unit.elementary_charge, sigma])
                    if net_charge > 0:
                        for w in r_waters:
                            w = int(w)
                            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(w)
                            wc = charge._value - 0.001
                            nonbonded_force.setParticleParameters(w, wc*unit.elementary_charge, sigma, epsilon)
                            for force in el_nb:
                                [l, charge, sigma] = force.getParticleParameters(w)
                                force.setParticleParameters(w, l, wc*unit.elementary_charge, sigma)



           #Calculate total energy for selected replicas

                if replica == i:

                    integrator_i = LangevinIntegrator(temperature, friction, dt)
                    integrator_i.setConstraintTolerance(constraintTolerance)
                    simulation_i = Simulation(topology, pH_system_temp, integrator_i, platform, platformProperties)
                    simulation_i.loadState('test-'+str(self._pH_list[i])+'.xml')
#                    simulation_i.minimizeEnergy(maxIterations = 1000)
                    positions_i = simulation_i.context.getState(getPositions=True).getPositions()

                elif replica == j:

                    integrator_j = LangevinIntegrator(temperature, friction, dt)
                    integrator_j.setConstraintTolerance(constraintTolerance)
                    simulation_j = Simulation(topology, pH_system_temp, integrator_j, platform, platformProperties)
                    simulation_j.loadState('test-'+str(self._pH_list[j])+'.xml')
#                    simulation_j.minimizeEnergy(maxIterations=1000)
                    positions_j = simulation_j.context.getState(getPositions=True).getPositions()
            simulation_i.context.setPositions(positions_j)
            simulation_j.context.setPositions(positions_i)
            state_ji = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            energy_ji = state_ji.getKineticEnergy()._value + state_ji.getPotentialEnergy()._value
            energy_ij = state_ij.getKineticEnergy()._value + state_ij.getPotentialEnergy()._value
            #Calculate total energy for replicas with exchnaged lambda values
            energy_ii_file = pd.read_csv('energy-'+str(self._pH_list[i])+'.csv')
            energy_ii = float(energy_ii_file.iloc[:, -1])
            energy_jj_file = pd.read_csv('energy-'+str(self._pH_list[j])+'.csv')
            energy_jj = float(energy_jj_file.iloc[:, -1])
            print("Replica exchange energy ", energy_ii, energy_jj, energy_ij, energy_ji)
            # Accept or reject the swap.
            log_p_accept = - (energy_ij + energy_ji) + energy_ii + energy_jj
            print(log_p_accept)
            if log_p_accept >= 0.0 or random.random() < math.exp(log_p_accept):
                print('Replica exchange accepted')
                # Swap states in replica slots i and j.
                simulation_i.context.setState(state_ij)
                simulation_i.reporters.append(CheckpointReporter('test-'+str(self._pH_list[i])+'.chk', 1))
                simulation_i.saveState('test-'+str(self._pH_list[i])+'.xml')
                simulation_i.context.createCheckpoint()
                simulation_j.context.setState(state_ji)
                simulation_j.reporters.append(CheckpointReporter('test-'+str(self._pH_list[j])+'.chk', 1))
                simulation_j.context.createCheckpoint()
                simulation_j.saveState('test-'+str(self._pH_list[j])+'.xml')
#                for residue in proton_change:
#                    proton = int(proton)
#                    residue = str(top.atom(proton)).split('-')[0]
#                    print('Shape ', lambda_list_i.shape[1], lambda_list_j.shape[1])
#                    lambda_i = liex.at[residue, str(lambda_list_i.shape[1])]
#                    lambda_j = ljex.at[residue, str(lambda_list_j.shape[1])]
#                    liex.at[residue, str(lambda_list_i.shape[1])] = lambda_j
#                    ljex.at[residue, str(lambda_list_j.shape[1])] = lambda_i
#                    df = pd.DataFrame.from_dict(lambda_list_i, orient='index', columns=[new_name])
#                    df_input = pd.read_csv('lambda_list-'+str(self._pH_list[i])+'.csv')
#                df_output = pd.concat([lambda_list_i, liex], axis = 1, sort = False)
#                df_output.to_csv('lambda_list-'+str(self._pH_list[i])+'.csv')
#                    df = pd.DataFrame.from_dict(lambda_list_j, orient='index', columns=[new_name])
#                    df_input = pd.read_csv('lambda_list-'+str(self._pH_list[j])+'.csv')
#                df_output = pd.concat([lambda_list_j, ljex], axis = 1, sort = False)
#                df_output.to_csv('lambda_list-'+str(self._pH_list[j])+'.csv')
            else:
                print('Replica_exchange_rejected')

parallel_tempering = pHRex(pH_system=pH_system, pH_list=pH_list)


# In[10]:


#parallel_tempering.run(n_iter = 1000, nsteps = 50000)
parallel_tempering.run(n_iter = 1, nsteps = 200)

