import sys, os
sys.path.append('../')

CWD = os.getcwd()
sys.path.append(CWD)

from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, calc_system_charge, create_force_particle, create_force_bond)
from assign_prot_prob import psf, list_alchem_residues, alchem_protons, segment_list
from input_file import *

########## intra_alchem_scaling = True is yet to be tested and added to the input_file. So far alchem-alchen non-bonded interactions are decoupled
intra_alchem_scaling = False
params = CharmmParameterSet(top_file, par_file)

CNB, CB = (np.array([], dtype=(np.int16)) for i in range(2))
system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, switchDistance=switchDistance)

pH_system = copy.deepcopy(system)


for force_index, reference_force in list(enumerate(pH_system.getForces())):
    reference_force_name = reference_force.__class__.__name__
    if reference_force_name == "NonbondedForce":
        NB = pH_system.getForces()[int(force_index)]

############################################################################
# FEP module (see OpenMM FEP)
#############################################################################

# Don't create a force if there are no alchemical atoms.
switch_width = 1 * angstroms

if len(list_alchem_residues) !=0:
    # --------------------------------------------------
    # Determine energy expression for all custom forces
    # --------------------------------------------------
    # Sterics mixing rules.
    sterics_mixing_rules = ('epsilon = sqrt(epsilon1*epsilon2);'  # Mixing rule for epsilon.
                            'sigma = 0.5*(sigma1 + sigma2);')  # Mixing rule for sigma.

    aa_sterics_mixing_rules = ('epsilon = sqrt(epsilon1*epsilon2);'
                               'sigma = 0.5*(sigma1 + sigma2);'
                               'lambda_sterics = lambda_sterics1*lambda_sterics2;')

    # Soft-core Lennard-Jones.
    exceptions_sterics_energy_expression = ('U_sterics;'
                                            'U_sterics = lambda_sterics*4*epsilon*x*(x-1.0);'
                                            'x = (sigma/reff_sterics)^6;'
                                             # Effective softcore distance for sterics.
                                            'reff_sterics = r;')


    # Define energy expression for sterics.
    sterics_energy_expression = exceptions_sterics_energy_expression + sterics_mixing_rules
    aa_sterics_energy_expression = exceptions_sterics_energy_expression + aa_sterics_mixing_rules

    # Define energy expression for electrostatics based on nonbonded method.
    nonbonded_method = NB.getNonbondedMethod()
    is_ewald_method = nonbonded_method in [openmm.NonbondedForce.Ewald,
                                           openmm.NonbondedForce.PME]
    is_rf_method = nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic,
                                        openmm.NonbondedForce.CutoffNonPeriodic]
    is_periodic_method = is_ewald_method or nonbonded_method == openmm.NonbondedForce.CutoffPeriodic

    energy_expressions = _get_electrostatics_energy_expressions(NB)
    (electrostatics_energy_expression, aa_electrostatics_energy_expression,
     exceptions_electrostatics_energy_expression) = energy_expressions  # Unpack tuple.

    # Create a copy of the NonbondedForce to handle particle interactions and
    # 1,4 exceptions between non-alchemical/non-alchemical atoms (nn).
    nonbonded_force = copy.deepcopy(NB)


    # Create CustomNonbondedForces to handle sterics particle interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). If intra-alchem_scaling is set to OFF, lambdas will be fixed
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_sterics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                     'lambda_sterics', is_lambda_controlled=True)
    if intra_alchem_scaling == True:
        aa_sterics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, aa_sterics_energy_expression,
                                                     'lambda_sterics', is_lambda_controlled=True)
    else:
        aa_sterics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                     'lambda_sterics', is_lambda_controlled=False)
    all_sterics_custom_nonbonded_forces = [na_sterics_custom_nonbonded_force, aa_sterics_custom_nonbonded_force]

    # Add parameters and configure CustomNonbondedForces to match reference force
    for force in all_sterics_custom_nonbonded_forces:
        force.addPerParticleParameter("lambda_sterics")
        force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
        force.addPerParticleParameter("epsilon")  # Lennard-Jones epsilon
        force.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction())
        force.setCutoffDistance(nonbonded_force.getCutoffDistance())
        force.setSwitchingDistance(nonbonded_force.getSwitchingDistance())
        force.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())
        if is_periodic_method:
            force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        else:
            force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

    # With exact PME treatment, particle electrostatics is handled through offset parameters.
    # Create CustomNonbondedForces to handle electrostatics particle interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_electrostatics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                            'lambda_electrostatics', is_lambda_controlled=True)
    if intra_alchem_scaling == True:
        aa_electrostatics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, aa_electrostatics_energy_expression,
                                                                'lambda_electrostatics', is_lambda_controlled=True)
    else:
        aa_electrostatics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                                'lambda_electrostatics', is_lambda_controlled=False)

    all_electrostatics_custom_nonbonded_forces = [na_electrostatics_custom_nonbonded_force,
                                                  aa_electrostatics_custom_nonbonded_force]

    # Common parameters and configuration for electrostatics CustomNonbondedForces.
    for force in all_electrostatics_custom_nonbonded_forces:
        force.addPerParticleParameter("lambda_electrostatics")
        force.addPerParticleParameter("charge")  # partial charge
        force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(nonbonded_force.getCutoffDistance() - switch_width)
        force.setCutoffDistance(nonbonded_force.getCutoffDistance())
        force.setUseLongRangeCorrection(False)  # long-range dispersion correction is meaningless for electrostatics

        if is_periodic_method:
            force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        else:
            force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

    # Create CustomBondForces to handle sterics 1,4 exceptions interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_sterics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                'lambda_sterics', is_lambda_controlled=True)
    aa_sterics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                'lambda_sterics', is_lambda_controlled=intra_alchem_scaling)
    all_sterics_custom_bond_forces = [na_sterics_custom_bond_force, aa_sterics_custom_bond_force]


    for force in all_sterics_custom_bond_forces:
        force.addPerBondParameter("lambda_sterics")
        force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma
        force.addPerBondParameter("epsilon")  # Lennard-Jones effective epsilon

    # Create CustomBondForces to handle electrostatics 1,4 exceptions interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_electrostatics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                       'lambda_electrostatics', is_lambda_controlled=True)
    aa_electrostatics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                       'lambda_electrostatics', is_lambda_controlled=intra_alchem_scaling)
    all_electrostatics_custom_bond_forces = [na_electrostatics_custom_bond_force, aa_electrostatics_custom_bond_force]

    # Create CustomBondForce to handle exceptions for electrostatics
    for force in all_electrostatics_custom_bond_forces:
        force.addPerBondParameter("lambda_electrostatics")
        force.addPerBondParameter("chargeprod")  # charge product
        force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma



   # -------------------------------------------------------------------------------
    # Distribute particle interactions contributions in appropriate nonbonded forces
    # -------------------------------------------------------------------------------

    # Create atom groups.
    alchemical_atomset = set(alchem_protons[0]) # all alchemical protons
    for segment in range(1, len(segment_list)):
        alchemical_atomset.update(set(alchem_protons[segment]))
    all_atomset = set(range(NB.getNumParticles()))  # all atoms, including alchemical region
    nonalchemical_atomset = all_atomset.difference(alchemical_atomset) #all non-alchemical atoms
    alchemical_atomset = all_atomset.difference(nonalchemical_atomset)

    # Fix any NonbondedForce issues with Lennard-Jones sigma = 0 (epsilon = 0), which should have sigma > 0.
    for particle_index in range(nonbonded_force.getNumParticles()):
        # Retrieve parameters.
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
        # Check particle sigma is not zero.
        if (sigma == 0.0 * angstrom):
            warning_msg = 'particle %d has Lennard-Jones sigma = 0 (charge=%s, sigma=%s, epsilon=%s); setting sigma=1A'
            logger.warning(warning_msg % (particle_index, str(charge), str(sigma), str(epsilon)))
            sigma = 1.0 * angstrom
            # Fix it.
            nonbonded_force.setParticleParameters(particle_index, charge, sigma, epsilon)

    for exception_index in range(nonbonded_force.getNumExceptions()):
        # Retrieve parameters.
        [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
        [charge_i, sigma, epsilon] = nonbonded_force.getParticleParameters(iatom)
        [charge_j, sigma, epsilon] = nonbonded_force.getParticleParameters(jatom)
        # Check particle sigma is not zero.
        if (sigma == 0.0 * angstrom):
            warning_msg = 'exception %d has Lennard-Jones sigma = 0 (iatom=%d, jatom=%d, chargeprod=%s, sigma=%s, epsilon=%s); setting sigma=1A'
            logger.warning(warning_msg % (exception_index, iatom, jatom, str(chargeprod), str(sigma), str(epsilon)))
            sigma = 1.0 * angstrom
            # Fix it.
            nonbonded_force.setExceptionParameters(exception_index, iatom, jatom, chargeprod, sigma, epsilon)
    # Copy NonbondedForce particle terms for alchemically-modified particles
    # to CustomNonbondedForces, and/or add the charge offsets for exact PME.
    # On CUDA, for efficiency reasons, all nonbonded forces (custom and not)
    # must have the same particles.
    for particle_index in range(nonbonded_force.getNumParticles()):
        # Retrieve nonbonded parameters.
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
        # Set sterics parameters in CustomNonbondedForces.
        for force in all_sterics_custom_nonbonded_forces:
            force.addParticle([1.0, sigma, epsilon])
            # Set electrostatics parameters in CustomNonbondedForces.
        for force in all_electrostatics_custom_nonbonded_forces:
            force.addParticle([1.0, charge, sigma])
    # NonbondedForce that will be handled by all other forces
    for particle_index in range(nonbonded_force.getNumParticles()):
        # Retrieve parameters.
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
        # Even with exact treatment of the PME electrostatics, we turn off
        # the NonbondedForce charge which is modeled by the offset parameter.
        if particle_index in alchemical_atomset:
            nonbonded_force.setParticleParameters(particle_index, abs(0.0*charge), sigma, abs(0*epsilon))
    # Restrict interaction evaluation of CustomNonbondedForces to their respective atom groups.
    na_sterics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomset)
    aa_sterics_custom_nonbonded_force.addInteractionGroup(alchemical_atomset, alchemical_atomset)
    na_electrostatics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomset)
    aa_electrostatics_custom_nonbonded_force.addInteractionGroup(alchemical_atomset, alchemical_atomset)

   # ---------------------------------------------------------------
    # Distribute exceptions contributions in appropriate bond forces
    # ---------------------------------------------------------------

    all_custom_nonbonded_forces = all_sterics_custom_nonbonded_forces + all_electrostatics_custom_nonbonded_forces

    # Move all NonbondedForce exception terms for alchemically-modified particles to CustomBondForces.
    for exception_index in range(nonbonded_force.getNumExceptions()):
        # Retrieve parameters.
        [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)

        # Exclude this atom pair in CustomNonbondedForces. All nonbonded forces
        # must have the same number of exceptions/exclusions on CUDA platform.
        for force in all_custom_nonbonded_forces:
            force.addExclusion(iatom, jatom)

        # Check how many alchemical atoms we have
        both_alchemical = iatom in alchemical_atomset and jatom in alchemical_atomset
        at_least_one_alchemical = iatom in alchemical_atomset or jatom in alchemical_atomset
        only_one_alchemical = at_least_one_alchemical and not both_alchemical

        # Check if this is an exception or an exclusion
        is_exception_epsilon = abs(epsilon.value_in_unit_system(md_unit_system)) > 0.0
        is_exception_chargeprod = abs(chargeprod.value_in_unit_system(md_unit_system)) > 0.0


        # If exception (and not exclusion), add special CustomBondForce terms to
        # handle alchemically-modified Lennard-Jones and electrostatics exceptions
        if both_alchemical:
            print('1-4 exception not handled for two alchemical atoms')

        if only_one_alchemical:
            if is_exception_epsilon:
                na_sterics_custom_bond_force.addBond(iatom, jatom, [1.0, sigma, epsilon])
            if is_exception_chargeprod:
                na_electrostatics_custom_bond_force.addBond(iatom, jatom, [1.0, chargeprod, sigma])

        # Turn off all exception contributions from alchemical atoms in the NonbondedForce
        # modelling non-alchemical atoms only
        if at_least_one_alchemical:
            nonbonded_force.setExceptionParameters(exception_index, iatom, jatom,
                                                   abs(0.0*chargeprod), sigma, abs(0.0*epsilon))

    # Add global parameters to forces.

    all_custom_forces = (all_custom_nonbonded_forces +
                         all_sterics_custom_bond_forces +
                         all_electrostatics_custom_bond_forces)


#Check if number of particles in all NB interections is the same

if na_electrostatics_custom_nonbonded_force.getNumParticles() == aa_electrostatics_custom_nonbonded_force.getNumParticles() == nonbonded_force.getNumParticles():
    for force_index, reference_force in list(enumerate(pH_system.getForces())):
        reference_force_name = reference_force.__class__.__name__
        if reference_force_name == "NonbondedForce":
            pH_system.removeForce(force_index)
            pH_system.addForce(nonbonded_force)
    for i in range (len(all_custom_forces)):
        pH_system.addForce(all_custom_forces[i])
else:
    print('Inconsistent number of nonbonded particles')

############################################################################
# Set a barostat according to input file
############################################################################

if barostat_type == "MonteCarlo":
    pH_system.addForce(MonteCarloBarostat(pressure, temperature))
elif barostat_type == "MonteCarloAnisotropic":
    pH_system.addForce(AnisotropicMonteCarloBarostat(pressure, temperature, scaleX = scale_X, scaleY = scale_Y, scaleZ = scale_Z, frequency = barostat_freq))
elif barostat_type == "MonteCarloMembrane":
    if xymode == 'XYIsotropic':
        if zmode == 'ZFree':
            pH_system.addForce(MonteCarloMembraneBarostat(pressure, surface_tension, temperature, MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
else:
    print("\nBarostat type not recognized\n")
    

for force_index, reference_force in list(enumerate(pH_system.getForces())):
    reference_force_name = reference_force.__class__.__name__
    if reference_force_name == "NonbondedForce":
        NBi = int(force_index)
    if reference_force_name == "CustomNonbondedForce" and reference_force.getNumPerParticleParameters() == 3:
        CNB = np.append(CNB, int(force_index))
    if reference_force_name == "CustomBondForce":
        CB = np.append(CB, int(force_index))
