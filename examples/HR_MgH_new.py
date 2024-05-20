import numpy as np
import sys
from helper_PFCI import PFHamiltonianGenerator
import scipy
import psi4
import json

mol_tmpl = """
1 1
Mg
H 1 **R**
symmetry c1
"""

mol_str = """
1 1
Mg
H 1 0.94
symmetry c1
"""

# number of bondlengths in scan
N_R = 80

# evenly spaced grid - can replace with a different grid (e.g. Chebyshev)
r_array = np.linspace(1, 3, N_R)


# number of CI roots to compute 
num_roots = 10  

# size of energy array will be num_roots x N_R
# energy_array = np.zeros((num_roots, N_R))

options_dict = {
        "basis": "6-31G",
        "scf_type": "pk",
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
}

lambda_list = [0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
for i in range(8):
        energy_array = np.zeros((num_roots, N_R))
        cavity_dict = {
        'omega_value' : 0.1941582620721647,
        'lambda_vector' : np.array([0, 0, lambda_list[i]]),
        'ci_level' : 'cas',
        'davidson_roots' : num_roots,
        'davidson_maxdim' : 13,
        'davidson_indim' : 3,
        'number_of_photons' : 5,
        'nact_orbs' : 12,
        'nact_els' : 8, 
        'photon_number_basis' : True,
        'canonical_mos' : True,
        'coherent_state_basis' : False
        }
        
        ctr = 0
        Dipoles = np.zeros((num_roots, num_roots, 3, N_R))
        for r in r_array:
                en_l = []
                mol_str = mol_tmpl.replace("**R**", str(r))
                mol = psi4.geometry(mol_str)
                print(mol_str)
                print(cavity_dict)
                print(options_dict)

                test_pf = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
                Dipoles[:,:,:,ctr] = np.copy(test_pf.dipole_array)
                energy_array[:,ctr] = np.copy(test_pf.CIeigs)
                ctr += 1

        np.save("fci_cavity_arrays_MgH_631G" + str(lambda_list[i]).replace(".", "_") , energy_array)
        np.save("dipoles_MgH"+str(lambda_list[i]).replace(".","_"), Dipoles)
        energy_lists = energy_array.tolist()
        r_list = r_array.tolist()
        dictionary = {
        "system": "MgH" + "_" +  str(lambda_list[i]).replace(".", "_"), # <== entered manually
        "r_data" : r_list, 
        "omega_value" : str(cavity_dict['omega_value']),  
        "basis_set" : str(options_dict["basis"]),
        "Photon basis" : "Number State", #<== Enter Manually
        "Number Photon States" : str(cavity_dict['number_of_photons']),
        "ci_level" : cavity_dict['ci_level'],
        "number of active orbitals" : str(cavity_dict['nact_orbs']), 
        "number of active electrons:" : str(cavity_dict['nact_els']),
        "energy_arrays" + str(lambda_list[i]).replace(".", "_") : energy_lists ,
        }
        json_object = json.dumps(dictionary, indent=4)

        with open(dictionary["system"] + ".json", "x") as outfile:
            outfile.write(json_object)

np.save("fci_r_array_MgH_6311G" + str(lambda_list[i]).replace(".", "_"), r_array)

