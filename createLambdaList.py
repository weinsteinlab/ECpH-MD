from imports import *
from fep_functions import (_get_pme_direct_space_unique_expression, _get_electrostatics_energy_expressions, add_global_parameters, calc_system_charge, create_force_particle, create_force_bond)
from definitions import *
from setup_pH_system import *

for pH in pH_list:
    lambda_list = {}
    lglu = float(l_Glu[pH])
    lhis = float(l_His[pH])
    llys = float(l_Lys[pH])
    lcys = float(l_Cys[pH])
    lasp = float(l_Asp[pH])
    for residue in list_alchem_residues:
        if l_special:
            for name in l_special:
                if residue == name:
                    lambda_list[residue] = l_special[name][pH]
        try: 
            print('residue with user-defined pKa ', residue, ' lambda value at pH ', pH, ': ', lambda_list[residue])
        except KeyError:
            if 'LYS' in residue:
                lambda_list[residue] = llys
            elif 'HSP' in residue:
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
