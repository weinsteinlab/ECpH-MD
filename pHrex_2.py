#!/usr/bin/env python
# coding: utf-8

# In[8]:


from imports import *
from definitions import *
from alchemical_protons import *
from fep_functions import *


def manage_waters(pH_system_temp):
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

def pick_charges(res_list, residue):
    for resname in res_list:
        if str(resname) in residue:
            charge_list = res_list[resname]
        elif 'HIS' in residue:
            if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                charge_list = res_list['HSE']
            else:
                charge_list = res_list['HSD']
    return charge_list

def create_cpH_system(pH_system_temp, lambda_list):
    nonbonded_force = pH_system_temp.getForces()[8]
    for side_atom in side_atoms:
        side_atom = int(side_atom)
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
        residue = str(top.atom(side_atom)).split('-')[0]
        atom_name = str(top.atom(side_atom)).split('-')[1]
        pick_charges(res_list, residue)
        for atom in charge_list:
            if atom_name in atom:
                n_ch = charge._value - charge_list[atom]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)
    for f in range(9, 17):
        force = pH_system_temp.getForces()[f]
        if force.__class__.__name__ == 'CustomNonbondedForce':
            for proton in (all_atoms):
                proton = int(proton)
                residue = str(top.atom(proton)).split('-')[0]
                [lambda_sterics, sigma, epsilon] = force.getParticleParameters(proton)
                if 'HIS' in residue:
                    atom_name = str(top.atom(proton)).split('-')[1]
                    if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 0.0:
                        if 'HE2' in atom_name:
                            lambda_sterics = lambda_list.at[residue, str(lambda_list.shape[1]-1)]
                        elif 'HD1' in atom_name:
                            lambda_sterics = 1.0
                    elif float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                        if 'HE2' in atom_name:
                            lambda_sterics = 1.0
                        elif 'HD1' in atom_name: 
                            lambda_sterics = lambda_list.at[residue, str(lambda_list.shape[1]-1)]
                else:        
                    lambda_sterics = lambda_list.at[residue, str(lambda_list.shape[1]-1)]

                force.setParticleParameters(proton, [lambda_sterics, sigma, epsilon])

            # Assign transition charge to neighboring alchemical atoms
            if force.getPerParticleParameterName(1) == 'charge':
                for side_atom in (side_atoms):
                    side_atom = int(side_atom)
                    [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                    residue = str(top.atom(side_atom)).split('-')[0]
                    atom_name = str(top.atom(side_atom)).split('-')[1]
                    pick_charges(res_list, residue)
                    for atom in charge_list:
                        if atom_name in atom:
                            n_ch = charge - charge_list[atom]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                             
        if force.__class__.__name__ == 'CustomBondForce':
            # Assign labmdas to alchemical protons                                               
            for bond in range(force.getNumBonds()):
                [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
                if atom_i in all_atoms:
                    proton = atom_i
                elif atom_j in all_atoms:
                    proton = atom_j
                residue = str(top.atom(proton)).split('-')[0]
                atom_name = str(top.atom(proton)).split('-')[1]
                if 'HIS' in residue:
                    if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 0.0:
                        if 'HE2' in atom_name:
                            a = lambda_list.at[residue, str(lambda_list.shape[1]-1)]
                        elif 'HD1' in atom_name:
                            a = 1.0
                    elif float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                        if 'HE2' in atom_name:
                            a = 1.0
                        elif 'HD1' in atom_name:
                            a = lambda_list.at[residue, str(lambda_list.shape[1]-1)]                            
                else:
                    a = lambda_list.at[residue, str(lambda_list.shape[1]-1)]

                force.setBondParameters(bond, q, r, [a, b, c])                                               

            # Assign transition charge to neighboring alchemical atoms
            if force.getPerBondParameterName(1) == 'chargeprod':
                for bond in range(force.getNumBonds()):
                    [atom_i, atom_j, (lambda_electrostatics, chargeprod, sigma)] = force.getBondParameters(bond)
                    if atom_i in side_atoms:
                        side_atom = atom_i
                        residue = str(top.atom(side_atom)).split('-')[0]
                        atom_name = str(top.atom(side_atom)).split('-')[1]
                        pick_charges(res_list, residue)
                        for atom in charge_list:
                            if atom_name in atom:
                                n_ch = b - charge_list[atom]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, n_ch, sigma]) 
                    elif atom_j in side_atoms:
                        side_atom = atom_j
                        residue = str(top.atom(side_atom)).split('-')[0]
                        atom_name = str(top.atom(side_atom)).split('-')[1]
                        pick_charges(res_list, residue)
                        for atom in charge_list:
                            if atom_name in atom:
                                n_ch = b - charge_list[atom]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, atom_i, atom_j, [lambda_electrostatics, n_ch, sigma])

                                
def create_cpH_exchange_system(pH_system, i, j, liex, ljex, new_name_i, new_name_j):
    pH_system_temp = copy.deepcopy(pH_system)
    rep_ex = [i, j]
    for replica in rep_ex:
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
            pick_charges(res_list, residue)
            for atom in charge_list:
                if atom_name in atom:
                    n_ch = charge._value - lambda_list[atom]*(1 - float(lambda_list.at[residue, new_name]))
                    nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)       
                    print(str(resname), ' parameters ', side_atom, n_ch, ' old charge: ', charge._value) 
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
                        pick_charges(res_list, residue)
                        for atom in charge_list:
                            if atom_name in atom:
                                n_ch = charge - charge_list[atom]*(1 - float(lambda_list.at[residue, new_name]))
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
                            pick_charges(res_list, residue)
                            for atom in charge_list:
                                if atom_name in atom:            
                                    n_ch = b - charge_list[atom]*(1 - float(lambda_list.at[residue, new_name]))
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

                            pick_charges(res_list, residue)
                            for atom in charge_list:
                                if atom_name in atom:
                                    n_ch = b - charge_list[atom]*(1 - float(lambda_list.at[residue, new_name]))
                                    force.setBondParameters(bond, q, r, [a, n_ch, c])  
        if replica == i:
            pH_system_temp_i = pH_system_temp
        elif replica == j:
            pH_system_temp_j = pH_system_temp
    return pH_system_temp_i, pH_system_temp_j



class pHrex:
    def __init__(self, pH_system, pH_list):
        self._pH_system = pH_system
        self._pH_list = pH_list
        
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
            #c = np.argwhere(self._pH_list==pH)[0][0]
            pH_system_temp = copy.deepcopy(self._pH_system)
            lambda_list = pd.read_csv('lambda_list-'+str(pH)+'.csv', index_col=0)
            create_cpH_system(pH_system_temp, lambda_list)
            manage_waters(pH_system_temp)
                                         
                                         
            integrator = LangevinIntegrator(temperature, friction, dt)
            integrator.setConstraintTolerance(constraintTolerance)
            simulation = Simulation(topology, pH_system_temp, integrator, platform, platformProperties)

            if iteration == 0:
                dcdReporter = DCDReporter('test-'+str(pH)+'.dcd', 1000)
                dataReporter = StateDataReporter('test-'+str(pH)+'.log', 100, totalSteps=steps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
                simulation.context.setPositions(positions_init)
                simulation.context.setVelocitiesToTemperature(temperature)
                simulation.currentStep = 0
                simulation.minimizeEnergy(maxIterations=1000)
                simulation.reporters.append(dcdReporter)
                simulation.reporters.append(dataReporter)
                simulation.step(nsteps)
                simulation.saveState('test-'+str(pH)+'.xml')
                state = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, 
                                                    getEnergy=True, getParameters=False, getParameterDerivatives=False)
                energy = state.getKineticEnergy()._value + state.getPotentialEnergy()._value
                energy = pd.DataFrame(np.array([energy]))
                energy.to_csv('energy-'+str(pH)+'.csv')
                simulation.minimizeEnergy(maxIterations = 1000)
                state_min = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False,
                                                    getEnergy=True, getParameters=False, getParameterDerivatives=False)
                energy_min = state_min.getKineticEnergy()._value + state_min.getPotentialEnergy()._value
                energy_min = pd.DataFrame(np.array([energy_min]))
                energy_min.to_csv('energy-min-'+str(pH)+'.csv')
            else:
                print('Iteration ', iteration, ' pH ', pH)  
                simulation.loadState('test-'+str(pH)+'.xml')           
                positions = simulation.context.getState(getPositions=True).getPositions()
                system_temp = simulation.context.getSystem()
             
                f = system_temp.getForces()[9]
                for i in all_atoms:
                    i = int(i)
                    print(i, ' parameters ', f.getParticleParameters(i))

                dcdReporter = DCDReporter('test-'+str(pH)+'.dcd', 1000, append = True)
                dataReporter = StateDataReporter('test-'+str(pH)+'.log', 100, totalSteps=steps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
                simulation.reporters.append(dcdReporter)
                simulation.reporters.append(dataReporter)
                simulation.step(nsteps)
                simulation.saveState('test-'+str(pH)+'.xml')
                state = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, 
                                                    getEnergy=True, getParameters=False, getParameterDerivatives=False)
                energy = state.getKineticEnergy()._value + state.getPotentialEnergy()._value
                energy_input = pd.read_csv('energy-'+str(pH)+'.csv', index_col = 0)
                energy_last = pd.DataFrame(np.array([energy]))
                name_i = str(energy_last.columns.tolist()[0])
                name_j = str(energy_input.shape[1])
                
                energy_current = energy_last.rename(columns={name_i:name_j})
                output = pd.concat([energy_input, energy_current], axis = 1, sort=False)
                output.to_csv('energy-'+str(pH)+'.csv')
                simulation.minimizeEnergy(maxIterations = 1000)
                state_min = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False,
                                                    getEnergy=True, getParameters=False, getParameterDerivatives=False)
                energy_min = state_min.getKineticEnergy()._value + state_min.getPotentialEnergy()._value
                energy_min_input = pd.read_csv('energy-min-'+str(pH)+'.csv', index_col = 0)
                energy_min_last = pd.DataFrame(np.array([energy_min]))
                name_i = str(energy_min_last.columns.tolist()[0])
                name_j = str(energy_min_input.shape[1])

                energy_min_current = energy_min_last.rename(columns={name_i:name_j})
                output = pd.concat([energy_min_input, energy_min_current], axis = 1, sort=False)
                output.to_csv('energy-min-'+str(pH)+'.csv')               



    def _mix_lambdas(self, n_attempts=1):
        # Attempt to switch two replicas at random. 
        for attempt in range(n_attempts):
            # Select two replicas at random.
            i = randint(0, self._pH_list.size-1)
            if i == 0:
                j = i+1
            else:
                j = i-1
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

            #pH_system_temp = copy.deepcopy(self._pH_system)            
            create_cpH_exchange_system(self._pH_system, i, j, liex, ljex, new_name_i, new_name_j)
            manage_waters(pH_system_temp)

            #Calculate total energy for selected replicas
            
            if replica == i:

                integrator_i = LangevinIntegrator(temperature, friction, dt)
                integrator_i.setConstraintTolerance(constraintTolerance)
                simulation_i = Simulation(topology, pH_system_temp_i, integrator_i, platform, platformProperties)
                simulation_i.loadState('test-'+str(self._pH_list[i])+'.xml')
                simulation_i.minimizeEnergy(maxIterations = 1000)

            elif replica == j:

                integrator_j = LangevinIntegrator(temperature, friction, dt)
                integrator_j.setConstraintTolerance(constraintTolerance)
                simulation_j = Simulation(topology, pH_system_temp_j, integrator_j, platform, platformProperties)
                simulation_j.loadState('test-'+str(self._pH_list[j])+'.xml')
                simulation_j.minimizeEnergy(maxIterations=1000)
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
                for residue in proton_change:
                    print('Shape ', lambda_list_i.shape[1], lambda_list_j.shape[1])
                    lambda_i = liex.at[residue, str(lambda_list_i.shape[1])]
                    lambda_j = ljex.at[residue, str(lambda_list_j.shape[1])]
                    liex.at[residue, str(lambda_list_i.shape[1])] = lambda_j
                    ljex.at[residue, str(lambda_list_j.shape[1])] = lambda_i
                df_output = pd.concat([lambda_list_i, liex], axis = 1, sort = False)
                df_output.to_csv('lambda_list-'+str(self._pH_list[i])+'.csv')
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
                if replica == i:
                    lambda_list = liex
                    new_name = new_name_i
                elif replica == j:
                    lambda_list = ljex
                    new_name = new_name_j
                else:
                    print("Error")
                create_cpH_system(pH_system_temp, lambda_list)    
                manage_waters(pH_system_temp)

           # Load the last coordinates and velocities of selected replicas

                if replica == i:

                    integrator_i = LangevinIntegrator(temperature, friction, dt)
                    integrator_i.setConstraintTolerance(constraintTolerance)
                    simulation_i = Simulation(topology, pH_system_temp, integrator_i, platform, platformProperties)
                    simulation_i.loadState('test-'+str(self._pH_list[i])+'.xml')
                    positions_i = simulation_i.context.getState(getPositions=True).getPositions()

                elif replica == j:

                    integrator_j = LangevinIntegrator(temperature, friction, dt)
                    integrator_j.setConstraintTolerance(constraintTolerance)
                    simulation_j = Simulation(topology, pH_system_temp, integrator_j, platform, platformProperties)
                    simulation_j.loadState('test-'+str(self._pH_list[j])+'.xml')
                    positions_j = simulation_j.context.getState(getPositions=True).getPositions()
                    
            # Exchange the posistions
            simulation_i.context.setPositions(positions_j)
            simulation_j.context.setPositions(positions_i)
            state_ji = simulation_j.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            state_ij = simulation_i.context.getState(getPositions=True, getVelocities=True, getForces=True,getEnergy=True, getParameters=True, getParameterDerivatives=True)
            
            # Calculate total energy for replicas with exchanged positions
            energy_ji = state_ji.getKineticEnergy()._value + state_ji.getPotentialEnergy()._value
            energy_ij = state_ij.getKineticEnergy()._value + state_ij.getPotentialEnergy()._value
            energy_ii_file = pd.read_csv('energy-'+str(self._pH_list[i])+'.csv')
            energy_ii = float(energy_ii_file.iloc[:, -1])
            energy_jj_file = pd.read_csv('energy-'+str(self._pH_list[j])+'.csv')
            energy_jj = float(energy_jj_file.iloc[:, -1])
            
            print("Replica exchange energy ", energy_ii, energy_jj, energy_ij, energy_ji)
            # Accept or reject the swap.
            log_p_accept = - (energy_ij + energy_ji) + energy_ii + energy_jj
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
            else:
                print('Replica_exchange_rejected')


# In[3]:





# In[ ]:




