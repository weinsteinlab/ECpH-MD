import sys,os
sys.path.append('../')

CWD = os.getcwd()
sys.path.append(CWD)


from imports import *
from input_file import *
ONE_4PI_EPS0 = 138.935456

def _get_pme_direct_space_unique_expression(reference_force):
    alpha_ewald, nx, ny, nz = reference_force.getPMEParameters()
    if alpha_ewald / alpha_ewald.unit == 0.0:
        tol = reference_force.getEwaldErrorTolerance()
        alpha_ewald = 1.0 / reference_force.getCutoffDistance() * np.sqrt(-np.log(2.0 * tol))
    alpha_ewald = alpha_ewald.value_in_unit_system(md_unit_system)
    pme_expression = '*erfc(alpha_ewald*reff_electrostatics)/reff_electrostatics;alpha_ewald = {};'.format(alpha_ewald)
    return pme_expression


def _get_electrostatics_energy_expressions(reference_force):
    electrostatics_prefix = 'U_electrostatics;U_electrostatics=lambda_electrostatics*ONE_4PI_EPS0*chargeprod'
    electrostatics_suffix = 'reff_electrostatics = r;ONE_4PI_EPS0 = {};'.format(ONE_4PI_EPS0)
    electrostatics_mixing_rules = 'chargeprod = charge1*charge2;sigma = 0.5*(sigma1 + sigma2);'
    aa_electrostatics_mixing_rules = 'chargeprod = charge1*charge2;sigma = 0.5*(sigma1 + sigma2); lambda_electrostatics = lambda_electrostatics1*lambda_electrostatics2'
    coulomb_expression = '/reff_electrostatics;'
    nonbonded_method = reference_force.getNonbondedMethod()
    if nonbonded_method in [openmm.NonbondedForce.NoCutoff]:
        electrostatics_method_expression = coulomb_expression
    elif nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic]:
        electrostatics_method_expression = self._get_reaction_field_unique_expression(reference_force)
    elif nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
        electrostatics_method_expression = _get_pme_direct_space_unique_expression(reference_force)
    else:
        raise ValueError('Nonbonded method {} not supported yet.'.format(nonbonded_method))
    exceptions_electrostatics_energy_expression = electrostatics_prefix
    exceptions_electrostatics_energy_expression += coulomb_expression
    exceptions_electrostatics_energy_expression += electrostatics_suffix
    electrostatics_energy_expression = electrostatics_prefix + electrostatics_method_expression + electrostatics_suffix + electrostatics_mixing_rules
    aa_electrostatics_energy_expression = electrostatics_prefix + electrostatics_method_expression + electrostatics_suffix + aa_electrostatics_mixing_rules
    return (
     electrostatics_energy_expression, aa_electrostatics_energy_expression, exceptions_electrostatics_energy_expression)



def calc_system_charge(force):
    net_charge = 0
    for i in range(force.getNumParticles()):
        lambdas, charge, epsilon = force.getParticleParameters(i)
        net_charge += round(charge, 3)

    net_charge = round(net_charge, 3)
    return net_charge


def create_force_particle(force_cls, energy_expression, lambda_variable_name, is_lambda_controlled):
    """Shortcut to create a lambda-controlled custom forces."""
    if is_lambda_controlled:
        force = force_cls(energy_expression)
        force.addGlobalParameter(lambda_variable_name, 1.0)
    else:
        energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
        force = force_cls(energy_expression_add)
    return force


def create_force_bond(force_cls, energy_expression, lambda_variable_name, is_lambda_controlled):
    """Shortcut to create a lambda-controlled custom forces."""
    if is_lambda_controlled:
        force = force_cls(energy_expression)
        force.addGlobalParameter(lambda_variable_name, 1.0)
    else:
        energy_expression_add = energy_expression + lambda_variable_name + '=1.0;'
        force = force_cls(energy_expression_add)
    return force
