import numpy as np
import sys
from helper_PFCI import PFHamiltonianGenerator
import scipy
import psi4

mol_tmpl = """
H
H 1 **R**
symmetry c1
"""

mol_str = """
H
H 1 0.74
symmetry c1
"""

# number of bondlengths in scan
N_R = 200

# evenly spaced grid - can replace with a different grid (e.g. Chebyshev)
r_array = np.linspace(0.5, 2.5, N_R)


# number of CI roots to compute 
num_roots = 12

# size of energy array will be num_roots x N_R
# energy_array = np.zeros((num_roots, N_R))

options_dict = {
        "basis": "6-311G",
        "scf_type": "pk",
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
}

lambda_list = [0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
for i in range(8):
        energy_array = np.zeros((num_roots, N_R))
        cavity_dict = {
        'omega_value' : 0.3944656448298627,
        'lambda_vector' : np.array([0, 0, lambda_list[i]]),
        'ci_level' : 'fci',
        'davidson_roots' : num_roots,
        'davidson_maxdim' : 13,
        'davidson_indim' : 12, 
        'number_of_photons' : 1,
        'photon_number_basis' : True,
        'canonical_mos' : True,
        'coherent_state_basis' : False
        }

        ctr = 0
        for r in r_array:
                en_l = []
                mol_str = mol_tmpl.replace("**R**", str(r))
                mol = psi4.geometry(mol_str)
                print(mol_str)
                print(cavity_dict)
                print(options_dict)

                test_pf = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
                energy_array[:,ctr] = np.copy(test_pf.CIeigs)
                ctr += 1

        np.save("fci_cavity_arrays_H2_6311G" + str(lambda_list[i]).replace(".", "_") , energy_array)
np.save("fci_r_array_H2_6311G" + str(lambda_list[i]).replace(".", "_"), r_array)




