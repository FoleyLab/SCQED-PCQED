import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator
np.set_printoptions(threshold=sys.maxsize)
psi4.core.set_output_file('output.dat', False)
import time
import json

qed_fci_tz_large_lambda = np.array([
 -2.9750971661e+00,
 -2.1822732987e+00,
 -2.0207189933e+00,
 -1.9983757383e+00,
 -1.7481245749e+00,
 -1.6950775096e+00,
 -1.6157398696e+00,
 -1.6157398696e+00,
 -1.5336748323e+00,
 -1.5336748323e+00]
)





# these file names should still be good
tz_en_file = "HHep_fci_cc_pVTZ_Energies.npy"
tz_mu_file = "HHep_fci_cc_pVTZ_Dipoles.npy"


tz_en = np.load(tz_en_file)
tz_mu = np.load(tz_mu_file)

# setup basic arguments to create an instance of the PFHamiltonianGenerator class
mol_str = """
Li
H 1 1.4
symmetry c1
"""

options_dict = {
    "basis": "sto-3g",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
}


cavity_free_dict = {
    'omega_value' : 0.0,
    'lambda_vector' : np.array([0, 0, 0.0]),
    'ci_level' : 'fci',   
    'full_diagonalization' : True,
    'number_of_photons' : 0, 
}

# create the instance of our PFHamiltonianGenerator class
tz_inst = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)

# number of adiabatic states for each basis set
N_el_tz = len(tz_en)

leng = 20

tz_dN = int(0.05 * N_el_tz)
print(F" Length of dz is {N_el_tz}, increment is {tz_dN}")


# array of number of electronic states in increments of 5% of total adiabatic states beginning from minimal basis
N_el_list_tz = np.linspace(tz_dN, leng * tz_dN, leng, dtype=int)
print(N_el_list_tz)


N_ph = 10
omega_dz = 0.9760568251
omega_tz = 0.9654959009
omega_qz = 0.9637811053



lambda_vector = np.array([0., 0., 0.02])

# double zeta arrays
tz_energies = np.zeros((leng, 10))
qed_fci_tz_energy = []

for i in range(leng):
    fast_start = time.time()
    # build double zeta Hamiltonian
    tz_inst.fast_build_pcqed_pf_hamiltonian(N_el_list_tz[i], N_ph, omega_tz, lambda_vector, tz_en, tz_mu, neglect_DSE=False)
    fast_end = time.time()
    dt = fast_end - fast_start
    ei = np.copy(tz_inst.PCQED_pf_eigs[:6])
    print(F' {N_el_list_tz[i]}, {ei[0]:12.10e},{ei[1]:12.10e},{ei[2]:12.10e},{ei[3]:12.10e},{ei[4]:12.10e},{ei[5]:12.10e}')
    # store double zeta energies
    tz_energies[i,:] = np.copy(tz_inst.PCQED_pf_eigs[:10])
    qed_fci_tz_energy.append(qed_fci_tz_large_lambda[0])


