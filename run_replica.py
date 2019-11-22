from imports import *
from input_file import *
from setup_pH_system import *
from pHrex import *
from definitions import *


pH = sys.argv[1]
iteration = int(sys.argv[2])
nsteps = int(sys.argv[3])

pH_system_temp = copy.deepcopy(pH_system)
lambda_list = pd.read_csv(('lambda_list-' + str(pH) + '.csv'), index_col=0)

create_cpH_system(pH_system_temp, lambda_list)
manage_waters(pH_system_temp)

# Integrator setup
if integrator_type == "Langevin":
    integrator = LangevinIntegrator(temperature, friction, dt)
elif integrator_type == "Verlet":
    integrator = VerletIntegrator(dt)
elif integrator_type == "Brownian":
    integrator = BrownianIntegrator(temperature, friction, dt)
elif integrator_type == "VariableVerlet":
    integrator = VariableVerletIntegrator(err_tol)
elif integrator_type == "VariableLangevin":
    integrator = VariableLangevinIntegrator(temperature, friction, err_tol)
else:
   print("\nIntegrator type not recognized!\n")

integrator.setConstraintTolerance(constraintTolerance)

# CUDA platform
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'0',  'Precision':'mixed'}

simulation = Simulation(topology, pH_system_temp, integrator, platform, platformProperties)

# Load state geometries 
if read_state_geometries:
    if len(pdb_state_files)== len(pH_list):
        i = np.where(pH_list == float(pH))[0][0]
        pdb = PDBFile(pdb_state_files[i])
        positions = pdb.positions
    else:
        print("\nThe pH-replicas and given states do not match\n")
else:
    pdb = PDBFile(pdb_file)
    positions = pdb.positions

if restart == 'OFF' and iteration == 0:
    print('Iteration ', iteration, 'pH ', pH)
    dcdReporter = DCDReporter(str(output_name) + '-' + str(pH) + '.dcd', dcdout_freq)
    dataReporter = StateDataReporter((str(output_name) + '-' + str(pH) + '.log'), 1000, totalSteps=nsteps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.currentStep = 0
    print('Minimizing...')
    simulation.minimizeEnergy(maxIterations=n_min_steps)
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    print('Simulating...')
    simulation.step(nsteps)
    simulation.saveState(str(output_name) + '-' + str(pH) + '-state.xml')
    state = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, getEnergy=True,
      getParameters=False,
      getParameterDerivatives=False)
    energy = state.getKineticEnergy()._value + state.getPotentialEnergy()._value
    energy = pd.DataFrame(np.array([energy]))
    energy.to_csv(str(output_name) + '-' + str(pH) + '-energy.csv')
    simulation.minimizeEnergy(maxIterations=1000)
    state_min = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, getEnergy=True,
      getParameters=False,
      getParameterDerivatives=False)
    energy_min = state_min.getKineticEnergy()._value + state_min.getPotentialEnergy()._value
    energy_min = pd.DataFrame(np.array([energy_min]))
    energy_min.to_csv(str(output_name) + '-' + str(pH) + '-energy-min.csv')
else:
    print('Iteration ', iteration, ' pH ', pH)
    simulation.loadState(str(output_name) + '-' + str(pH) + '-state.xml')
    positions = simulation.context.getState(getPositions=True).getPositions()
    system_temp = simulation.context.getSystem()
    dcdReporter = DCDReporter((str(output_name) + '-' + str(pH) + '.dcd'), dcdout_freq, append=True)
    dataReporter = StateDataReporter((str(output_name) + '-' + str(pH) + '.log'), 1000, totalSteps=nsteps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    print('Simulating...')
    simulation.step(nsteps)
    simulation.saveState(str(output_name) + '-' + str(pH) + '-state.xml')
    state = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, getEnergy=True,
      getParameters=False,
      getParameterDerivatives=False)
    energy = state.getKineticEnergy()._value + state.getPotentialEnergy()._value
    energy_input = pd.read_csv((str(output_name) + '-' + str(pH) + '-energy.csv'), index_col=0)
    energy_last = pd.DataFrame(np.array([energy]))
    name_i = str(energy_last.columns.tolist()[0])
    name_j = str(energy_input.shape[1])
    energy_current = energy_last.rename(columns={name_i: name_j})
    output = pd.concat([energy_input, energy_current], axis=1, sort=False)
    output.to_csv(str(output_name) + '-' + str(pH) + '-energy.csv')
    simulation.minimizeEnergy(maxIterations=1000)
    state_min = simulation.context.getState(getPositions=False, getVelocities=False, getForces=False, getEnergy=True,
      getParameters=False,
      getParameterDerivatives=False)
    energy_min = state_min.getKineticEnergy()._value + state_min.getPotentialEnergy()._value
    energy_min_input = pd.read_csv((str(output_name) + '-' + str(pH) + '-energy-min.csv'), index_col=0)
    energy_min_last = pd.DataFrame(np.array([energy_min]))
    name_i = str(energy_min_last.columns.tolist()[0])
    name_j = str(energy_min_input.shape[1])
    energy_min_current = energy_min_last.rename(columns={name_i: name_j})
    output = pd.concat([energy_min_input, energy_min_current], axis=1, sort=False)
    output.to_csv(str(output_name) + '-' + str(pH) + '-energy-min.csv')
