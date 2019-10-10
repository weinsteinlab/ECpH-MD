#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import glob
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout
#from openmmtools import alchemy
#from openmmtools.alchemy import *
#from openmmtools.alchemy import AbsoluteAlchemicalFactory
import math
from math import pi
from openmmtools import testsystems, alchemy
from simtk import openmm, unit
from openmmtools.states import *
#from openmmtools import testsystems, states, mcmc
import simtk.openmm.app.dcdfile
#from openmmtools import testsystems, alchemy
from simtk import openmm, unit
import pandas as pd
import mdtraj as md
import numpy as np
import random
from random import randint
from openmmtools import cache
from openmmtools.constants import ONE_4PI_EPS0


# In[2]:


def _get_pme_direct_space_unique_expression(reference_force):
    # Determine PME parameters.
    [alpha_ewald, nx, ny, nz] = reference_force.getPMEParameters()
    if (alpha_ewald/alpha_ewald.unit) == 0.0:
        # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance.
        tol = reference_force.getEwaldErrorTolerance()
        alpha_ewald = (1.0/reference_force.getCutoffDistance()) * np.sqrt(-np.log(2.0*tol))
    alpha_ewald = alpha_ewald.value_in_unit_system(unit.md_unit_system)
    pme_expression = ("*erfc(alpha_ewald*reff_electrostatics)/reff_electrostatics;"
                      "alpha_ewald = {};").format(alpha_ewald)
    return pme_expression

def _get_electrostatics_energy_expressions(reference_force):
    # The final expression will be prefix + method + suffix.
    electrostatics_prefix = ('U_electrostatics;'
                             'U_electrostatics=(lambda_electrostatics^softcore_d)*ONE_4PI_EPS0*chargeprod')
    # Effective softcore distance for electrostatics (common to all methods).
    electrostatics_suffix = ('reff_electrostatics = sigma*((softcore_beta*(1.0-lambda_electrostatics)^softcore_e + (r/sigma)^softcore_f))^(1/softcore_f);'
                             'ONE_4PI_EPS0 = {};').format(ONE_4PI_EPS0)  # Already in OpenMM units.
    # Define mixing rules.
    electrostatics_mixing_rules = ('chargeprod = charge1*charge2;'  # Mixing rule for charges.
                                   'sigma = 0.5*(sigma1 + sigma2);')  # Mixing rule for sigma.
    # Standard Coulomb expression with softened core. This is used
    coulomb_expression = '/reff_electrostatics;'
    # Select electrostatics functional form based on nonbonded method.
    nonbonded_method = reference_force.getNonbondedMethod()
    # Soft-core Coulomb.
    if nonbonded_method in [openmm.NonbondedForce.NoCutoff]:
        electrostatics_method_expression = coulomb_expression
    # Reaction-field electrostatics.
    elif nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic]:
        electrostatics_method_expression = self._get_reaction_field_unique_expression(reference_force)
    # PME electrostatics.
    elif nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
        # Ewald direct-space electrostatics.
        electrostatics_method_expression = _get_pme_direct_space_unique_expression(reference_force)
    else:
        raise ValueError("Nonbonded method {} not supported yet.".format(nonbonded_method))
        # Define energy expression for 1,4 electrostatic exceptions.
    exceptions_electrostatics_energy_expression = electrostatics_prefix
     #    exceptions_electrostatics_energy_expression += electrostatics_method_expression
    exceptions_electrostatics_energy_expression += coulomb_expression
    exceptions_electrostatics_energy_expression += electrostatics_suffix
    # Define energy expression for electrostatics.
    electrostatics_energy_expression = (electrostatics_prefix + electrostatics_method_expression +
                                        electrostatics_suffix + electrostatics_mixing_rules)
    return electrostatics_energy_expression, exceptions_electrostatics_energy_expression

def add_global_parameters(force):
    force.addGlobalParameter('softcore_alpha', softcore_alpha)
    force.addGlobalParameter('softcore_beta', softcore_beta)
    force.addGlobalParameter('softcore_a', softcore_a)
    force.addGlobalParameter('softcore_b', softcore_b)
    force.addGlobalParameter('softcore_c', softcore_c)
    force.addGlobalParameter('softcore_d', softcore_d)
    force.addGlobalParameter('softcore_e', softcore_e)
    force.addGlobalParameter('softcore_f', softcore_f)
def calc_system_charge(force):
    c = 0
    for i in range(force.getNumParticles()):
        [laas, charge, epsilon] = force.getParticleParameters(i)
        c+=round(charge, 3)
    c = round(c, 3)
    print('Net charge: ', c)
    return c


# In[3]:


print('Building system...')
psf = CharmmPsfFile('../fep-40-i01-1-wb-i.psf')
pdb = PDBFile('../fep-45-i01-1-wb-i.pdb')
topology = psf.topology
positions_init = pdb.positions
params = CharmmParameterSet('../all_top.rtf', '../parameters.prm')
nonbondedMethod = PME
nonbondedCutoff = 12*angstroms
switchDistance=10*angstroms

ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001

# Integration Options

dt = 0.002*picoseconds
temperature = 276*kelvin
friction = 1.0/picosecond
#integrator = LangevinIntegrator(temperature, friction, dt)
#integrator.setConstraintTolerance(constraintTolerance)

x_PBC_vector_length = 6.05
y_PBC_vector_length = 5.88
z_PBC_vector_length = 5.25
  
psf.setBox(x_PBC_vector_length, y_PBC_vector_length, z_PBC_vector_length)

system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, switchDistance=switchDistance)
pH_system = copy.deepcopy(system)

platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex': '0', 'Precision': 'mixed'}


# In[8]:


t = md.load('../fep-45-i01-1-wb-i.pdb')
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
#asp_side_atoms = top.select('rescode ASP2 and name CB CG OD1 OD2')
side_atoms = np.concatenate((lys_side_atoms, his_side_atoms, glu_side_atoms))
side_atoms = np.sort(side_atoms)

# Define VdW softcore parameters
switch_width=1*unit.angstroms
softcore_alpha = 0.5
softcore_a = 1
softcore_b = 1
softcore_c = 1
softcore_beta = 0.0
softcore_d = 1
softcore_e = 1
softcore_f = 1
# Define the replica exchange procedure
pH_low = 2
pH_high = 12
pH_step = 0.5
steps = 100
# Default pKa values
EpKa = 4.4
KpKa = 10.4
HpKa = 6.5
HpKa_switch = 9.1
CpKa = 9.5
DpKa = 4.0
#Lambda lists
l_Glu = {}
l_His = {}
l_Lys = {}
l_Cys = {}
l_Asp = {}
l_Glu_special = {}
glu_char = {}
lys_char = {}
his_char_E = {}
his_char_D = {}
asp_char = {}
cys_char = {}
# GLU difference charges
glu_char["CG:"] = 0.07
glu_char["CD:"] = 0.13
glu_char["OE1:"] = 0.21
glu_char["OE2:"] = 0.15
# LYS difference charges
lys_char["CE:"] = 0.08
lys_char["HE1:"] = -0.025
lys_char["HE2:"] = -0.025
lys_char["NZ:"] = 0.66
lys_char["HZ1:"] = -0.01
lys_char["HZ2:"] = -0.01
# HSE difference charges
his_char_E["CB:"] = 0.03
his_char_E["ND1:"] = 0.19
#his_char_E["HD1:"] = 0.44
his_char_E["CG:"] = -0.03
his_char_E["CE1:"] = 0.07
his_char_E["HE1:"] = 0.05
his_char_E["NE2:"] = -0.15
his_char_E["HE2:"] = 0.12
his_char_E["CD2:"] = 0.24
his_char_E["HD2:"] = 0.04
# HSD difference charges
his_char_D["CB:"] = 0.04
his_char_D["ND1:"] = -0.15
his_char_D["HD1:"] = 0.12
his_char_D["CG:"] = 0.24
his_char_D["CE1:"] = 0.07
his_char_D["HE1:"] = 0.05
his_char_D["NE2:"] = 0.19
#his_char_D["HE2:"] = 0.44
his_char_D["CD2:"] = -0.03
his_char_D["HD2:"] = 0.03
# CYS difference charges
cys_char["CB:"] = 0.27
cys_char["SG:"] = 0.57
# ASP difference charges
asp_char["CB:"] = 0.07
asp_char["CG:"] = 0.13
asp_char["OD1:"] = 0.21
asp_char["OD2:"] = 0.15

for i in range(pH_low*10, pH_high*10, 5):
    p = 1-1/(1+10**(EpKa-i/10))
    l_Glu[i/10] = '%.4f'%p
    p = 1-1/(1+10**(HpKa-i/10))
    l_His[i/10] = '%.4f'%p
    p = 1-1/(1+10**(KpKa-i/10))
    l_Lys[i/10] = '%.4f'%p
    p = 1-1/(1+10**(CpKa-i/10))
    l_Cys[i/10] = '%.4f'%p
    p = 1-1/(1+10**(DpKa-i/10))
    l_Asp[i/10] = '%.4f'%p
pi = round(pi, 5)

pH_list = np.arange(pH_low, pH_high, 0.5)
print(pH_list)
all_atoms_HSE = np.concatenate((lys_atoms, his_atoms_E, glu_atoms))
all_atoms_HSE = np.sort(all_atoms_HSE)
all_atoms = np.concatenate((lys_atoms, his_atoms_E, his_atoms_D, glu_atoms))
all_atoms = np.sort(all_atoms)
list_alchem_residues = [None]*all_atoms_HSE.size

for i in range(len(list_alchem_residues)):
    list_alchem_residues[i] = str(top.atom(all_atoms_HSE[i])).split('-')[0]

for pH in pH_list:
    lambda_list = {}
    lglu = float(l_Glu[pH])
    lhis = float(l_His[pH])
    llys = float(l_Lys[pH])
    lcys = float(l_Cys[pH])
    lasp = float(l_Asp[pH])
    for residue in list_alchem_residues:
        if 'LYS' in residue:
            lambda_list[residue] = llys
        elif 'HIS' in residue:
            lambda_list[residue] = lhis
            if pH < HpKa_switch:
                lambda_list[residue+'_sw'] = 1
            else:
                lambda_list[residue+'_sw'] = 0
        elif 'GLU' in residue:
            lambda_list[residue] = lglu
        elif 'ASP' in residue:
            lambda_list[residue] = lasp
        elif 'CYS' in residue:
            lambda_list[residue] = lcys
    df = pd.DataFrame.from_dict(lambda_list, orient='index')
    df.to_csv('lambda_list-'+str(pH)+'.csv')
print(lambda_list)


# In[9]:


for force_index, reference_force in list(enumerate(pH_system.getForces())):
    reference_force_name = reference_force.__class__.__name__
    if reference_force_name == "NonbondedForce":
        NB = pH_system.getForces()[int(force_index)]
    elif reference_force_name == "AndersenThermostat":
        pH_system.removeForce(int(force_index))
#print('pH = ', pH, 'lambdas: GLU ', lglu, ' LYS ', llys, ' HIS ', lhis, ' CYS ', lcys, ' ASP ', lasp)


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
    def create_force_particle(force_cls, energy_expression, lambda_variable_name, is_lambda_controlled):
        """Shortcut to create a lambda-controlled custom forces."""
        if is_lambda_controlled:
            energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
            force = force_cls(energy_expression_add)
            #force.addPerParticleParameter(lambda_variable_name)
        else:  # fix lambda variable to 1.0
            energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
            force = force_cls(energy_expression_add)
        return force
    def create_force_bond(force_cls, energy_expression, lambda_variable_name, is_lambda_controlled):
        """Shortcut to create a lambda-controlled custom forces."""
        if is_lambda_controlled:
            energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
            force = force_cls(energy_expression_add)

        else:  # fix lambda variable to 1.0
            energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
            force = force_cls(energy_expression_add)
        return force

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
    #alchemical_atomset_lys = lys_atoms
    #alchemical_atomset_glu = glu_atoms
    #alchemical_atomset_his = his_atoms
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
            #if particle_index in lys_atoms:
            #    force.addParticle([llys, sigma, epsilon])
            #elif particle_index in his_atoms:
            #    force.addParticle([lhis, sigma, epsilon])
            #elif particle_index in glu_atoms:
            #    force.addParticle([lglu, sigma, epsilon])
            #elif particle_index in cys_atoms:
            #    force.addParticle([lcys, sigma, epsilon])
            #elif particle_index in asp_atoms:
            #    force.addParticle([lasp, sigma, epsilon])                
            #else:
            force.addParticle([1.0, sigma, epsilon])
            # Set electrostatics parameters in CustomNonbondedForces.
        for force in all_electrostatics_custom_nonbonded_forces:
            #if particle_index in lys_atoms:
            #    force.addParticle([llys, charge, sigma])
            #elif particle_index in his_atoms:
            #    force.addParticle([lhis, charge, sigma])
            #elif particle_index in glu_atoms:
            #    force.addParticle([lglu, charge, sigma]) 
            #elif particle_index in cys_atoms:
            #    force.addParticle([lcys, charge, sigma])
            #elif particle_index in asp_atoms:
            #    force.addParticle([lasp, charge, sigma])                
            #else: 
            force.addParticle([1.0, charge, sigma])
                #if 'LYS'in str(psf.atom_list[particle_index]):
                #    for i in lys_char:
                #        if i in str(psf.atom_list[particle_index]):
                #            #print(str(psf.atom_list[particle_index]))
                #            npch = charge._value - (1-llys)*lys_char[i]
                #            force.setParticleParameters(particle_index, [1.0, npch, sigma])                           
                #elif 'HSP'in str(psf.atom_list[particle_index]):
                #    for i in his_char:
                #        if i in str(psf.atom_list[particle_index]):
                #            #print(str(psf.atom_list[particle_index]))
                #            npch = charge._value - (1-lhis)*his_char[i]
                #            force.setParticleParameters(particle_index, [1.0, npch, sigma])
                #elif 'GLU'in str(psf.atom_list[particle_index]):
                #    for i in glu_char:
                #        if i in str(psf.atom_list[particle_index]):
                #            #print(str(psf.atom_list[particle_index]))
                #            npch = charge._value - (1 - lglu)*glu_char[i]
                #            force.setParticleParameters(particle_index, [1.0, npch, sigma])     
                #elif 'CYS'in str(psf.atom_list[particle_index]):
                #    for i in cys_char:
                #        if i in str(psf.atom_list[particle_index]):
                #            #print(str(psf.atom_list[particle_index]))
                #            npch = charge._value - (1 - lcys)*cys_char[i]
                #            force.setParticleParameters(particle_index, [1.0, npch, sigma]) 
                #elif 'ASP2'in str(psf.atom_list[particle_index]):
                #    for i in asp_char:
                #        if i in str(psf.atom_list[particle_index]):
                #            #print(str(psf.atom_list[particle_index])
                #            npch = charge._value - (1-lasp)*asp_char[i]
                #            force.setParticleParameters(particle_index, [1.0, npch, sigma])                             
        # Turn off interactions contribution from alchemically-modified particles in unmodified
    # NonbondedForce that will be handled by all other forces
    for particle_index in range(nonbonded_force.getNumParticles()):
        # Retrieve parameters.
        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
        # Even with exact treatment of the PME electrostatics, we turn off
        # the NonbondedForce charge which is modeled by the offset parameter.
        if particle_index in alchemical_atomset:
            nonbonded_force.setParticleParameters(particle_index, abs(0.0*charge), sigma, abs(0*epsilon))
        #if 'LYS'in str(psf.atom_list[particle_index]):
        #    for i in lys_char:
        #        if i in str(psf.atom_list[particle_index]):
        #            npch = charge._value - (1-llys)*lys_char[i]
        #            nonbonded_force.setParticleParameters(particle_index, npch, sigma, epsilon)
        #if 'HSP'in str(psf.atom_list[particle_index]):
        #    for i in his_char:
        #        if i in str(psf.atom_list[particle_index]):
        #            npch = charge._value - (1-lhis)*his_char[i]
        #            nonbonded_force.setParticleParameters(particle_index, npch, sigma, epsilon)                    
        #if 'GLU2'in str(psf.atom_list[particle_index]):
        #    for i in glu_char:
        #        if i in str(psf.atom_list[particle_index]):
        #            npch = charge._value - (1-lglu)*glu_char[i]
        #            nonbonded_force.setParticleParameters(particle_index, npch, sigma, epsilon)
        #if 'CYS'in str(psf.atom_list[particle_index]):
        #    for i in cys_char:
        #        if i in str(psf.atom_list[particle_index]):
        #            npch = charge._value - (1-lcys)*cys_char[i]
        #            nonbonded_force.setParticleParameters(particle_index, npch, sigma, epsilon)
        #if 'ASP2'in str(psf.atom_list[particle_index]):
        #    for i in asp_char:
        #        if i in str(psf.atom_list[particle_index]):
        #            npch = charge._value - (1-lasp)*asp_char[i]
        #            nonbonded_force.setParticleParameters(particle_index, npch, sigma, epsilon)                    
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
                #if iatom or jatom in lys_atoms:    
                #    na_sterics_custom_bond_force.addBond(iatom, jatom, [llys, sigma, epsilon])
                #elif iatom or jatom in his_atoms:
                #    na_sterics_custom_bond_force.addBond(iatom, jatom, [lhis, sigma, epsilon])
                #elif iatom or jatom in glu_atoms:
                #    na_sterics_custom_bond_force.addBond(iatom, jatom, [lglu, sigma, epsilon])
                #elif iatom or jatom in cys_atoms:
                #    na_sterics_custom_bond_force.addBond(iatom, jatom, [lcys, sigma, epsilon])
                #elif iatom or jatom in asp_atoms:
                na_sterics_custom_bond_force.addBond(iatom, jatom, [1.0, sigma, epsilon])                     
            if is_exception_chargeprod:
                #if iatom or jatom in lys_atoms:    
                #    na_electrostatics_custom_bond_force.addBond(iatom, jatom, [llys, chargeprod, sigma])
                #elif iatom or jatom in his_atoms:
                #    na_electrostatics_custom_bond_force.addBond(iatom, jatom, [lhis, chargeprod, sigma])
                #elif iatom or jatom in glu_atoms:
                #    na_electrostatics_custom_bond_force.addBond(iatom, jatom, [lglu, chargeprod, sigma])
                #elif iatom or jatom in cys_atoms:
                #    na_electrostatics_custom_bond_force.addBond(iatom, jatom, [lcys, chargeprod, sigma])
                #elif iatom or jatom in asp_atoms:
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


# In[22]:


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
            c = np.argwhere(self._pH_list==pH)[0][0]
            pH_system_temp = copy.deepcopy(self._pH_system)
            lambda_list = pd.read_csv('lambda_list-'+str(pH)+'.csv', index_col=0)
            nonbonded_force = pH_system_temp.getForces()[8]
            for side_atom in side_atoms:
                side_atom = int(side_atom)
                [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(side_atom)
                residue = str(top.atom(side_atom)).split('-')[0]
                atom_name = str(top.atom(side_atom)).split('-')[1]
                if 'LYS' in residue:
                    for charge_delta in lys_char:
                        if atom_name in charge_delta:
                            n_ch = charge._value - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                            print(charge, ' ', n_ch)
                            nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)               
                if 'GLU' in residue:
                    for charge_delta in glu_char:
                        if atom_name in charge_delta:
                            n_ch = charge._value - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                            print(charge, ' ', n_ch)
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
#                            print(charge, ' ', n_ch)
                            nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)      
                if 'CYS' in residue:
                    for charge_delta in cys_char:
                        if atom_name in charge_delta:
                            n_ch = charge._value - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                            print(charge, ' ', n_ch)
                            nonbonded_force.setParticleParameters(side_atom, n_ch, sigma, epsilon)       
            for f in range(9, 17):
                force = pH_system_temp.getForces()[f]
                if force.__class__.__name__ == 'CustomNonbondedForce':
                    # Assign labmdas to alchemical protons
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

                    # Assign transition charge to neighboring alchemical atoms
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
#                                        print(charge, ' ', n_ch)
                                        force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma])                        
                            elif 'GLU' in residue:
                                for charge_delta in glu_char:
                                    if atom_name in charge_delta:
                                        n_ch = charge - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                        print(charge, ' ', n_ch)
                                        force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                            elif 'HIS' in residue:
                                if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                                    his_char = his_char_E
                                else:
                                    his_char = his_char_D
                                for charge_delta in his_char:
                                    if atom_name in charge_delta:
                                        n_ch = charge - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                        print(charge, ' ', n_ch)
                                        force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                            elif 'ASP' in residue:
                                for charge_delta in asp_char:
                                    if atom_name in charge_delta:
                                        n_ch = charge - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                        print(charge, ' ', n_ch)
                                        force.setParticleParameters(side_atom, [lambda_electrostatics, n_ch, sigma]) 
                            elif 'CYS' in residue:
                                for charge_delta in cys_char:
                                    if atom_name in charge_delta:
                                        n_ch = charge - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
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
                            [q, r, (a, b, c)] = force.getBondParameters(bond)
                            #b = round(b, 5)
                            if q in side_atoms:
                                side_atom = q
                                residue = str(top.atom(side_atom)).split('-')[0]
                                atom_name = str(top.atom(side_atom)).split('-')[1]
                                if 'LYS' in residue:
                                    for charge_delta in lys_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])                                              
                                if 'GLU' in residue:
                                    for charge_delta in glu_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'HIS' in residue:
                                    if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                                        his_char = his_char_E
                                    else:
                                        his_char = his_char_D
                                    for charge_delta in his_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'ASP' in residue:
                                    for charge_delta in asp_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'CYS' in residue:
                                    for charge_delta in cys_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c]) 
                            elif r in side_atoms:
                                side_atom = r
                                residue = str(top.atom(side_atom)).split('-')[0]
                                atom_name = str(top.atom(side_atom)).split('-')[1]
                                if 'LYS' in residue:
                                    for charge_delta in lys_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - lys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'GLU' in residue:
                                    for charge_delta in glu_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - glu_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'HIS' in residue:
                                    if float(lambda_list.at[residue+'_sw', str(lambda_list.shape[1]-1)]) == 1.0:
                                        his_char = his_char_E
                                    else:
                                        his_char = his_char_D
                                    for charge_delta in his_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - his_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'ASP' in residue:
                                    for charge_delta in asp_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - asp_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
#                                            print(b, ' ', n_ch)
                                            force.setBondParameters(bond, q, r, [a, n_ch, c])
                                if 'CYS' in residue:
                                    for charge_delta in cys_char:
                                        if atom_name in charge_delta:
                                            n_ch = b - cys_char[charge_delta]*(1 - float(lambda_list.at[residue, str(lambda_list.shape[1]-1)]))
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
                            force.setParticleParameters(i, l, wc*unit.elementary_charge, sigma)               
            integrator = LangevinIntegrator(temperature, friction, dt)
            integrator.setConstraintTolerance(constraintTolerance)
            simulation = Simulation(topology, pH_system_temp, integrator, platform, platformProperties)
#            if abs(net_charge) > 0.0001:
#                n_waters = int(abs(net_charge)/0.001)
#                print('Number of modified water oxygens ', n_waters)
#
#                r_waters = list(np.random.choice(waters, n_waters))
#                if net_charge < 0:
#                    for i in r_waters:
#                        i = int(i)            
#                        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(i)
#                        wc = charge._value + 0.001
#                        nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
#                        for force in el_nb:
#                            [l, charge, sigma] = force.getParticleParameters(i)
#                            force.setParticleParameters(i, [l, wc*unit.elementary_charge, sigma])
#                if net_charge > 0:
#                    for i in r_waters:
#                        i = int(i)            
#                        [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(i)
#                        wc = charge._value - 0.001        
#                        nonbonded_force.setParticleParameters(i, wc*unit.elementary_charge, sigma, epsilon)
#                        for force in el_nb:
#                            [l, charge, sigma] = force.getParticleParameters(i)
#                            force.setParticleParameters(i, l, wc*unit.elementary_charge, sigma)                    
#            simulation.context.reinitialize()


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
#                energy_last = energy_input.iloc[:, -1]
#                energy_last = pd.DataFrame(energy_last)
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
#                energy_last = energy_input.iloc[:, -1]
#                energy_last = pd.DataFrame(energy_last)
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


parallel_tempering.run(n_iter = 1000, nsteps = 50000)

