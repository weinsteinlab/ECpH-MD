from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
from definitions import *

from setup_pH_system import *
from pHrex_2 import *

pH = sys.argv[1]
iteration = int(sys.argv[2])
nsteps = int(sys.argv[3])

pH_system_temp = pH_system
lambda_list = pd.read_csv('lambda_list-'+str(pH)+'.csv', index_col=0)
create_cpH_system(pH_system_temp, lambda_list)
manage_waters(pH_system_temp)


integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, pH_system_temp, integrator, platform, platformProperties)
#simulation = Simulation(topology, pH_system_temp, integrator)

if iteration == 0:
    print('akjsfkajghskfjghk')
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
    print('oiujsfhkmnxbv,asjhd')
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

