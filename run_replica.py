from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
from definitions import *

from setup_pH_system import *

pH = sys.argv[1]
iteration = int(sys.argv[2])
nsteps = int(sys.argv[3])

pH_system_temp = pH_system
lambda_list = pd.read_csv('lambda_list-'+str(pH)+'.csv', index_col=0)
nonbonded_force = pH_system_temp.getForces()[8]

# Assign the transition charges to atoms neighboring to alchemical protons according to current pH value for non-alchemical interactions

for side_atom in side_atoms:
    side_atom = int(side_atom)
    [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
    residue = str(top.atom(side_atom)).split('-')[0]
    atom_name = str(top.atom(side_atom)).split('-')[1]
    if 'LYS' in residue:
        for charge_delta in lys_char:
            if atom_name in charge_delta:
                n_ch = charge._value - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)               
    if 'GLU' in residue:
        for charge_delta in glu_char:
            if atom_name in charge_delta:
                n_ch = charge._value - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
    if 'HIS' in residue:
        if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
            his_char = his_char_E
        else:
            his_char = his_char_D
        for charge_delta in his_char:
            if atom_name in charge_delta:
                n_ch = charge._value - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
    if 'ASP' in residue:
        for charge_delta in asp_char:
            if atom_name in charge_delta:
                n_ch = charge._value - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
    if 'CYS' in residue:
        for charge_delta in cys_char:
            if atom_name in charge_delta:
                n_ch = charge._value - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)       
for f in range(9, 17):
    force = pH_system_temp.getForces()[f]
    if force.__class__.__name__ == 'CustomNonbondedForce':

        # Assign lambdas to alchemical protons

        for proton in (all_atoms):
            proton = int(proton)
            residue = str(top.atom(proton)).split('-')[0]
            [a, b, c] = force.getParticleParameters(proton)
            if 'HIS' in residue:
                atom_name = str(top.atom(proton)).split('-')[1]
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

            force.setParticleParameters(proton, [a, b, c])

        # Assign the transition charges to atoms neighboring to alchemical protons according to current pH value for alchemical interactions

        if force.getPerParticleParameterName(1) == 'charge':
            for side_atom in (side_atoms):
                side_atom = int(side_atom)
                [lambda_electrostatics, charge, sigma] = force.getParticleParameters(side_atom)
                residue = str(top.atom(side_atom)).split('-')[0]
                atom_name = str(top.atom(side_atom)).split('-')[1]
                if 'LYS' in residue:
                    for charge_delta in lys_char:
                        if atom_name in charge_delta:
                            n_ch = charge - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                        
                elif 'GLU' in residue:
                    for charge_delta in glu_char:
                        if atom_name in charge_delta:
                            n_ch = charge - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                elif 'HIS' in residue:
                    if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                        his_char = his_char_E
                    else:
                        his_char = his_char_D
                    for charge_delta in his_char:
                        if atom_name in charge_delta:
                            n_ch = charge - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                elif 'ASP' in residue:
                    for charge_delta in asp_char:
                        if atom_name in charge_delta:
                            n_ch = charge - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                elif 'CYS' in residue:
                    for charge_delta in cys_char:
                        if atom_name in charge_delta:
                            n_ch = charge - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                            force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                             
    if force.__class__.__name__ == 'CustomBondForce':

        # Assign lambdas to alchemical protons 

        for bond in range(force.getNumBonds()):
            [q, r, (a, b, c)] = force.getBondParameters(bond)
            if q in all_atoms:
                proton = q
            elif r in all_atoms:
                proton = r
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

        # Assign the transition charges to atoms neighboring to alchemical protons according to current pH value for alchemical interactions

        if force.getPerBondParameterName(1) == 'chargeprod':
            for bond in range(force.getNumBonds()):
                [q, r, (a, b, c)] = force.getBondParameters(bond)
                if q in side_atoms:
                    side_atom = q
                    residue = str(top.atom(side_atom)).split('-')[0]
                    atom_name = str(top.atom(side_atom)).split('-')[1]
                    if 'LYS' in residue:
                        for charge_delta in lys_char:
                            if atom_name in charge_delta:
                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])                                              
                    if 'GLU' in residue:
                        for charge_delta in glu_char:
                            if atom_name in charge_delta:
                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'HIS' in residue:
                        if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                            his_char = his_char_E
                        else:
                            his_char = his_char_D
                        for charge_delta in his_char:
                            if atom_name in charge_delta:
                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'ASP' in residue:
                        for charge_delta in asp_char:
                            if atom_name in charge_delta:
                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'CYS' in residue:
                        for charge_delta in cys_char:
                            if atom_name in charge_delta:
                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c]) 
                elif r in side_atoms:
                    side_atom = r
                    residue = str(top.atom(side_atom)).split('-')[0]
                    atom_name = str(top.atom(side_atom)).split('-')[1]
                    if 'LYS' in residue:
                        for charge_delta in lys_char:
                            if atom_name in charge_delta:
                                n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'GLU' in residue:
                        for charge_delta in glu_char:
                            if atom_name in charge_delta:
                                n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'HIS' in residue:
                        if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                            his_char = his_char_E
                        else:
                            his_char = his_char_D
                        for charge_delta in his_char:
                            if atom_name in charge_delta:
                                n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'ASP' in residue:
                        for charge_delta in asp_char:
                            if atom_name in charge_delta:
                                n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])
                    if 'CYS' in residue:
                        for charge_delta in cys_char:
                            if atom_name in charge_delta:
                                n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
                                force.setBondParameters(bond, q, r, [a, n_ch, c])

# Calculate the net charge of the system

el_nb = [pH_system_temp.getForces()[11], pH_system_temp.getForces()[12]]
net_charge = calc_system_charge(el_nb[0])
# If the net chatge is higher than a given value alter the partial charges of random water oxygens by 0.001 e
if abs(net_charge) > 0.001:
    n_waters = int(abs(net_charge)/0.001)
    print('Number of modified water oxygens ', n_waters)
    # if the charge of the water oxygen has been changed previously restore it 
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
        for i in r_waters:
            i = int(i)            
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(i)
            wc = charge._value + 0.001
            nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
            for force in el_nb:
                [l, charge, sigma] = force.getParticleParameters(i)
                force.setParticleParameters(i, [l, wc*unit.elementary_charge, sigma])
    if net_charge > 0:
        for i in r_waters:
            i = int(i)            
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(i)
            wc = charge._value - 0.001        
            nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
            for force in el_nb:
                [l, charge, sigma] = force.getParticleParameters(i)
                force.setParticleParameters(i, l, wc*unit.elementary_charge, siigma)               

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
 
    # For dubugging mode print all the parameters of all alchemical protons
   
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
