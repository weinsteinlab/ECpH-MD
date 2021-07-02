import sys, os
sys.path.append('../')

CWD = os.getcwd()
sys.path.append(CWD)

from imports import *
from assign_prot_prob import  pH_list, l_Glu, l_His, l_Lys, l_Cys, l_Asp, list_alchem_residues, l_special, accuracy
from definitions import HpKa_switch
for pH in np.round(pH_list, accuracy):
    lambda_list = {}
    lglu = float(l_Glu[pH])
    lhis = float(l_His[pH])
    llys = float(l_Lys[pH])
    lcys = float(l_Cys[pH])
    lasp = float(l_Asp[pH])
    for residue in list_alchem_residues:
        if l_special:
            for name in l_special:
                if name in residue:
                    lambda_list[residue] = l_special[name][pH]
        try: 
            print('\nresidue with user-defined pKa: ', residue, '\nlambda value at pH ', pH, ': ', lambda_list[residue])
        except KeyError:
            if 'LYS' in residue:
                lambda_list[residue] = llys
            elif 'HSP' in residue or 'HIS' in residue:
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
    df.to_csv('./lambdas/lambda_list-'+str(pH)+'.csv')
