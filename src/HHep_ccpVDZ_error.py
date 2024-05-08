import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator
np.set_printoptions(threshold=sys.maxsize)
psi4.core.set_output_file('output.dat', False)
import time
import json

qed_fci_dz_large_lambda = np.array([ 
 -2.9607270006e+00,
 -2.1651918650e+00,
 -1.9953443012e+00,
 -1.9738708169e+00,
 -1.6473292910e+00,
 -1.5828421987e+00,
 -1.1892368461e+00,
 -1.1235525952e+00,
 -1.1235525952e+00,
 -1.0954722116e+00]
)





# these file names should still be good
dz_en_file = "HHep_fci_cc_pVDZ_Energies.npy"
dz_mu_file = "HHep_fci_cc_pVDZ_Dipoles.npy"


dz_en = np.load(dz_en_file)
dz_mu = np.load(dz_mu_file)

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
dz_inst = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)

# number of adiabatic states for each basis set
N_el_dz = len(dz_en)

beg = 10
leng = 20

dz_dN = int(0.05 * N_el_dz)
print(F" Length of dz is {N_el_dz}, increment is {dz_dN}")


# array of number of electronic states in increments of 5% of total adiabatic states beginning from minimal basis
N_el_list_dz = np.linspace(dz_dN, leng * dz_dN, leng, dtype=int)
print(N_el_list_dz)


N_ph = 10
omega_dz = 0.9760568251
omega_tz = 0.9654959009
omega_qz = 0.9637811053



lambda_vector = np.array([0., 0., 0.02])

# double zeta arrays
dz_energies = np.zeros((leng, 10))
qed_fci_dz_energy = []

for i in range(beg,leng):
    fast_start = time.time()
    # build double zeta Hamiltonian
    dz_inst.fast_build_pcqed_pf_hamiltonian(N_el_list_dz[i], N_ph, omega_dz, lambda_vector, dz_en, dz_mu, neglect_DSE=False)
    fast_end = time.time()
    dt = fast_end - fast_start
    ei = np.copy(dz_inst.PCQED_pf_eigs[:6])
    print(F' {N_el_list_dz[i]}, {ei[0]:12.10e},{ei[1]:12.10e},{ei[2]:12.10e},{ei[3]:12.10e},{ei[4]:12.10e},{ei[5]:12.10e}')
    # store double zeta energies
    dz_energies[i,:] = np.copy(dz_inst.PCQED_pf_eigs[:10])
    qed_fci_dz_energy.append(qed_fci_dz_large_lambda[0])


