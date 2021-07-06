import sys, os

CWD = os.getcwd()
sys.path.append(CWD)

from imports import *
from input_file import *
from setup_pH_system import *
from pHrex import *
#from definitions import *

pH = sys.argv[1]
subjob_number = int(sys.argv[2])
replica_number = int(sys.argv[3])

print("============================================================================================================")
print('======pH:' + str(pH) + ' subjob_number:' + str(subjob_number) + ' replica_number:' + str(replica_number)) 
print("============================================================================================================")
print()

lambda_list = pd.read_csv(('./lambdas/lambda_list-' + str(pH) + '.csv'), index_col=0)
create_cpH_system(pH_system, lambda_list)
manage_waters(pH_system)

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
platformProperties = {'DeviceIndex':'0',  'Precision':'mixed'}

simulation = Simulation(topology, pH_system, integrator, platform, platformProperties)

# Load structure 

pdb = PDBFile(pdb_file)
positions = pdb.positions


print(' pH:', pH, ' subjob:', subjob_number, 'replica:', replica_number)

output_directory = './simulations/pH_' + str(pH) + '_replica_number_' + str(replica_number).zfill(4) + '/'
output_no_path = str(output_name) + '-ph' + str(pH) + '_replica_number_' + str(replica_number).zfill(4) + '-subjob' + str(subjob_number).zfill(4) 
full_output_name = output_directory + str(output_no_path)

if subjob_number == 0:
    print("Reading in previous state: FALSE")
    dcdReporter = DCDReporter(full_output_name + '.dcd', dcdout_freq)
    dataReporter = StateDataReporter((full_output_name + '.log'), 10000, totalSteps=md_steps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.currentStep = 0

    print('Minimizing...')
    simulation.minimizeEnergy(maxIterations=n_min_steps)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(full_output_name+'-min.pdb', 'w'))
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)

    print('Simulating...')
    simulation.step(md_steps)
    print('Simulation completed\nSaving current state...')
    simulation.saveState(full_output_name + '-state.xml')

    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
    ### In current version of the code with OpenMM 7.5.0 the units of potential energy (kJ/mol) and Boltzmann constant (J/mol do not match
    energy = state.getPotentialEnergy()/(AVOGADRO_CONSTANT_NA*BOLTZMANN_CONSTANT_kB*temperature*1000)
    energy = pd.DataFrame(np.array([energy]))
    energy.to_csv('./energies/' + str(output_no_path) + '-energy.csv')
    print('FINISH')

else:
    previous_output_name=output_directory + str(output_name) + '-ph' + str(pH) + '_replica_number_' + str(replica_number).zfill(4) + '-subjob' + str(subjob_number - 1).zfill(4)
    print("Reading previous state from: " + previous_output_name + "-state.xml")
    simulation.loadState(previous_output_name + '-state.xml')
    positions = simulation.context.getState(getPositions=True).getPositions()
    system_temp = simulation.context.getSystem()
    dcdReporter = DCDReporter((full_output_name + '.dcd'), dcdout_freq)
    dataReporter = StateDataReporter((full_output_name + '.log'), 10000, totalSteps=md_steps, step=True, time=True, speed=True, progress=True, elapsedTime=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, separator=',')
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)

    print('Simulating...')
    simulation.step(md_steps)
    print('Simulation completed\nSaving current state...')
    simulation.saveState(full_output_name + '-state.xml')

    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, getParameterDerivatives=True)
    ### In current version of the code with OpenMM 7.5.0 the units of potential energy (kJ/mol) and Boltzmann constant (J/mol do not match
    energy = state.getPotentialEnergy()/(AVOGADRO_CONSTANT_NA*BOLTZMANN_CONSTANT_kB*temperature*1000) 
    energy_input = pd.read_csv(('./energies/' + str(output_name) + '-ph' + str(pH) + '_replica_number_' + str(replica_number).zfill(4) +     '-subjob' + str(subjob_number - 1).zfill(4) + '-energy.csv'), index_col=0)
    energy_last = pd.DataFrame(np.array([energy]))
    name_i = str(energy_last.columns.tolist()[0])
    name_j = str(energy_input.shape[1])
    energy_current = energy_last.rename(columns={name_i: name_j})
    output = pd.concat([energy_input, energy_current], axis=1, sort=False)
    output.to_csv('./energies/' + str(output_no_path) + '-energy.csv')
    print('FINISH')
