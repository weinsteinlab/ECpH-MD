from imports import *
from pHrex import *
from pH_exchange_functions import *
from definitions import *
from buildSystem import *

system = buildSystem()
pH_system = copy.deepcopy(system)

t = md.load('./inputFiles/fep-45-i01-1-wb-i.pdb')
top = t.topology

# Define alchemical protons
lys_atoms = top.select('resname LYS and name HZ3')
his_atoms_E = top.select('rescode H and name HD1')
his_atoms_D = top.select('rescode H and name HE2')
glu_atoms = top.select('resname GLU2 and name HE2')
asp_atoms = top.select('resname ASP2 and name HD2')
cys_atoms = top.select('resname CYS and name HG1')

# Select water molecules to keep the net charge
waters = t.top.select('water and type O')

#Define side atoms with flexible charge
lys_side_atoms = top.select('rescode K and name CE  HE2 HE3 NZ HZ1 HZ2')
his_side_atoms = top.select('rescode H and name CB ND1 CG CE1 HE1 NE2 CD2 HD2 HE2 HD1')
glu_side_atoms = top.select('resname GLU2 and name CG CD OE1 OE2')

side_atoms = np.concatenate((lys_side_atoms, his_side_atoms, glu_side_atoms))
side_atoms = np.sort(side_atoms)


print(pH_list)
all_atoms_HSE = np.concatenate((lys_atoms, his_atoms_E, glu_atoms))
all_atoms_HSE = np.sort(all_atoms_HSE)
all_atoms = np.concatenate((lys_atoms, his_atoms_E, his_atoms_D, glu_atoms))
all_atoms = np.sort(all_atoms)
list_alchem_residues = [None]*all_atoms_HSE.size

for i in range(len(list_alchem_residues)):
    list_alchem_residues[i] = str(top.atom(all_atoms_HSE[i])).split('-')[0]


for force_index, reference_force in list(enumerate(pH_system.getForces())):
    reference_force_name = reference_force.__class__.__name__
    if reference_force_name == "NonbondedForce":
        NB = pH_system.getForces()[int(force_index)]
    elif reference_force_name == "AndersenThermostat":
        pH_system.removeForce(int(force_index))


# Don't create a force if there are no alchemical atoms.
if len(list_alchem_residues) !=0:
    # --------------------------------------------------
    # Determine energy expression for all custom forces
    # --------------------------------------------------
    # Sterics mixing rules.
    sterics_mixing_rules = ('epsilon = sqrt(epsilon1*epsilon2);'  # Mixing rule for epsilon.
                            'sigma = 0.5*(sigma1 + sigma2);')  # Mixing rule for sigma.

    # Soft-core Lennard-Jones.
    exceptions_sterics_energy_expression = ('U_sterics;'
                                            'U_sterics = (lambda_sterics^softcore_a)*4*epsilon*x*(x-1.0);'
                                            'x = (sigma/reff_sterics)^6;'
                                             # Effective softcore distance for sterics.
                                            'reff_sterics = sigma*((softcore_alpha*(1.0-lambda_sterics)^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);')

    # Define energy expression for sterics.
    sterics_energy_expression = exceptions_sterics_energy_expression + sterics_mixing_rules

    # Define energy expression for electrostatics based on nonbonded method.
    nonbonded_method = NB.getNonbondedMethod()
    is_ewald_method = nonbonded_method in [openmm.NonbondedForce.Ewald,
                                           openmm.NonbondedForce.PME]
    is_rf_method = nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic,
                                        openmm.NonbondedForce.CutoffNonPeriodic]
    is_periodic_method = is_ewald_method or nonbonded_method == openmm.NonbondedForce.CutoffPeriodic

    energy_expressions = _get_electrostatics_energy_expressions(NB)
    (electrostatics_energy_expression,
     exceptions_electrostatics_energy_expression) = energy_expressions  # Unpack tuple.
    # ------------------------------------------------------------
    # Create and configure all forces to add to alchemical system
    # ------------------------------------------------------------
    # Interactions and exceptions will be distributed according to the following table.
    # --------------------------------------------------------------------------------------------------
    # FORCE                                    | INTERACTION GROUP                                     |
    # --------------------------------------------------------------------------------------------------
    # nonbonded_force (unmodified)             | all interactions nonalchemical/nonalchemical          |
    #                                          | all exceptions nonalchemical/nonalchemical            |
    # --------------------------------------------------------------------------------------------------
    # aa_sterics_custom_nonbonded_force        | sterics interactions alchemical/alchemical            |
    # --------------------------------------------------------------------------------------------------
    # aa_electrostatics_custom_nonbonded_force | electrostatics interactions alchemical/alchemical     |
    #                                          | (only without exact PME treatment)                    |
    # --------------------------------------------------------------------------------------------------
    # na_sterics_custom_nonbonded_force        | sterics interactions non-alchemical/alchemical        |
    # --------------------------------------------------------------------------------------------------
    # na_electrostatics_custom_nonbonded_force | electrostatics interactions non-alchemical/alchemical |
    #                                          | (only without exact PME treatment)                    |
    # --------------------------------------------------------------------------------------------------
    # aa_sterics_custom_bond_force             | sterics exceptions alchemical/alchemical              |
    # --------------------------------------------------------------------------------------------------
    # aa_electrostatics_custom_bond_force      | electrostatics exceptions alchemical/alchemical       |
    #                                          | (only without exact PME treatment)                    |
    # --------------------------------------------------------------------------------------------------
    # na_sterics_custom_bond_force             | sterics exceptions non-alchemical/alchemical          |
    # --------------------------------------------------------------------------------------------------
    # na_electrostatics_custom_bond_force      | electrostatics exceptions non-alchemical/alchemical   |
    #                                          | (only without exact PME treatment)                    |
    # --------------------------------------------------------------------------------------------------

    # Create a copy of the NonbondedForce to handle particle interactions and
    # 1,4 exceptions between non-alchemical/non-alchemical atoms (nn).
    nonbonded_force = copy.deepcopy(NB)

    # Create CustomNonbondedForces to handle sterics particle interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_sterics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, sterics_energy_expression,
                                                     'lambda_sterics', is_lambda_controlled=True)
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
    aa_electrostatics_custom_nonbonded_force = create_force_particle(openmm.CustomNonbondedForce, electrostatics_energy_expression,
                                                            'lambda_electrostatics', is_lambda_controlled=False)
    all_electrostatics_custom_nonbonded_forces = [na_electrostatics_custom_nonbonded_force,
                                                  aa_electrostatics_custom_nonbonded_force]

    # Common parameters and configuration for electrostatics CustomNonbondedForces.
    for force_el_alchem in all_electrostatics_custom_nonbonded_forces:
        force_el_alchem.addPerParticleParameter("lambda_electrostatics")
        force_el_alchem.addPerParticleParameter("charge")  # partial charge
        force_el_alchem.addPerParticleParameter("sigma")  # Lennard-Jones sigma
        force_el_alchem.setUseSwitchingFunction(True)
        force_el_alchem.setSwitchingDistance(nonbonded_force.getCutoffDistance() - switch_width)
        force_el_alchem.setCutoffDistance(nonbonded_force.getCutoffDistance())
        force_el_alchem.setUseLongRangeCorrection(False)  # long-range dispersion correction is meaningless for electrostatics

        if is_periodic_method:
            force_el_alchem.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        else:
            force_el_alchem.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

    # Create CustomBondForces to handle sterics 1,4 exceptions interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_sterics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                'lambda_sterics', is_lambda_controlled=True)
    aa_sterics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_sterics_energy_expression,
                                                'lambda_sterics', is_lambda_controlled=False)
    all_sterics_custom_bond_forces = [na_sterics_custom_bond_force, aa_sterics_custom_bond_force]

    for force_blj_alchem in all_sterics_custom_bond_forces:
        force_blj_alchem.addPerBondParameter("lambda_sterics")
        force_blj_alchem.addPerBondParameter("sigma")  # Lennard-Jones effective sigma
        force_blj_alchem.addPerBondParameter("epsilon")  # Lennard-Jones effective epsilon

    # Create CustomBondForces to handle electrostatics 1,4 exceptions interactions between
    # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
    # to 1.0 for decoupled interactions in alchemical/alchemical force.
    na_electrostatics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                       'lambda_electrostatics', is_lambda_controlled=True)
    aa_electrostatics_custom_bond_force = create_force_bond(openmm.CustomBondForce, exceptions_electrostatics_energy_expression,
                                                       'lambda_electrostatics', is_lambda_controlled=False)
    all_electrostatics_custom_bond_forces = [na_electrostatics_custom_bond_force, aa_electrostatics_custom_bond_force]

    # Create CustomBondForce to handle exceptions for electrostatics
    for force_bel_alchem in all_electrostatics_custom_bond_forces:
        force_bel_alchem.addPerBondParameter("lambda_electrostatics")
        force_bel_alchem.addPerBondParameter("chargeprod")  # charge product
        force_bel_alchem.addPerBondParameter("sigma")  # Lennard-Jones effective sigma

    # -------------------------------------------------------------------------------
    # Distribute particle interactions contributions in appropriate nonbonded forces
    # -------------------------------------------------------------------------------
    # Create atom groups.
    alchemical_atomset = set(all_atoms)
    all_atomset = set(range(NB.getNumParticles()))  # all atoms, including alchemical region
    nonalchemical_atomset = all_atomset.difference(alchemical_atomset)
    alchemical_atomset = all_atomset.difference(nonalchemical_atomset)


    # Fix any NonbondedForce issues with Lennard-Jones sigma = 0 (epsilon = 0), which should have sigma > 0.
    for particle_index in range(nonbonded_force.getNumParticles()):
        # Retrieve parameters.
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
        # Check particle sigma is not zero.
        if (sigma == 0.0 * unit.angstrom):
            warning_msg = 'particle %d has Lennard-Jones sigma = 0 (charge=%s, sigma=%s, epsilon=%s); setting sigma=1A'
            logger.warning(warning_msg % (particle_index, str(charge), str(sigma), str(epsilon)))
            sigma = 1.0 * unit.angstrom
            # Fix it.
            nonbonded_force.setParticleParameters(particle_index, charge, sigma, epsilon)


    for exception_index in range(nonbonded_force.getNumExceptions()):
        # Retrieve parameters.
        [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
        [charge_i, sigma, epsilon] = nonbonded_force.getParticleParameters(iatom)                    
        [charge_j, sigma, epsilon] = nonbonded_force.getParticleParameters(jatom)                                  
        # Check particle sigma is not zero.
        if (sigma == 0.0 * unit.angstrom):
            warning_msg = 'exception %d has Lennard-Jones sigma = 0 (iatom=%d, jatom=%d, chargeprod=%s, sigma=%s, epsilon=%s); setting sigma=1A'
            logger.warning(warning_msg % (exception_index, iatom, jatom, str(chargeprod), str(sigma), str(epsilon)))
            sigma = 1.0 * unit.angstrom
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
        for force in all_electrostatics_custom_nonbonded_forces:
            force.addParticle([1.0, charge, sigma])
        # Turn off interactions contribution from alchemically-modified particles in unmodified


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
        is_exception_epsilon = abs(epsilon.value_in_unit_system(unit.md_unit_system)) > 0.0
        is_exception_chargeprod = abs(chargeprod.value_in_unit_system(unit.md_unit_system)) > 0.0


        # If exception (and not exclusion), add special CustomBondForce terms to
        # handle alchemically-modified Lennard-Jones and electrostatics exceptions
        if both_alchemical:
            print('1-4 exception not handled for two alchemical atoms')

        if only_one_alchemical:
            if is_exception_epsilon:
                na_sterics_custom_bond_force.addBond(iatom, jatom, [1.0, sigma, epsilon])                     
            if is_exception_chargeprod:
                na_electrostatics_custom_bond_force.addBond(iatom, jatom, [1.0, chargeprod, sigma])                    
        # else: both particles are non-alchemical, leave them in the unmodified NonbondedForce

        # Turn off all exception contributions from alchemical atoms in the NonbondedForce
        # modelling non-alchemical atoms only
        if at_least_one_alchemical:
            nonbonded_force.setExceptionParameters(exception_index, iatom, jatom,
                                                   abs(0.0*chargeprod), sigma, abs(0.0*epsilon))

    # Add global parameters to forces.
    all_custom_forces = (all_custom_nonbonded_forces +
                         all_sterics_custom_bond_forces +
                         all_electrostatics_custom_bond_forces)
    for force in all_custom_forces:
        add_global_parameters(force)


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


pH_system.addForce(MonteCarloBarostat(1.01325*bar, 276*kelvin))
pH_system.getForces()
parallel_tempering = pHrex(pH_system=pH_system, pH_list=pH_list)
parallel_tempering.run(n_iter = 1000, nsteps = 1000)
