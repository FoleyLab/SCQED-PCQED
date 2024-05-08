import numpy as np
import sys
from helper_PFCI import PFHamiltonianGenerator
import scipy
import psi4

mol_tmpl = """
Li
C 1 **R**
N 2 1.0 1 180
symmetry c1
"""

mol_str = """
Li
C 1 1.0
N 2 1.0 1 180
symmetry c1
"""

# number of bondlengths in scan
N_R = 500

# evenly spaced grid - can replace with a different grid (e.g. Chebyshev)
r_array = np.linspace(0.5, 3.0, N_R)


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
        'omega_value' : 0.2460433494247667,
        'lambda_vector' : np.array([0, 0, lambda_list[i]]),
        'ci_level' : 'cas',
        'davidson_roots' : num_roots,
        'davidson_maxdim' : 13,
        'davidson_indim' : 3,
        'nact_orbs' : 13,
        'nact_els' : 14, 
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

        np.save("cis_cavity_arrays_LiCN_sto3g" + str(lambda_list[i]).replace(".", "_") , energy_array)
np.save("cis_r_array_LiCN", r_array)




