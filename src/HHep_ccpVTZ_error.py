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

qed_fci_qz_large_lambda = np.array([
 -2.9773960939e+00,
 -2.1846350116e+00,
 -2.0247781044e+00,
 -2.0023579319e+00,
 -1.7700528337e+00,
 -1.7270874847e+00,
 -1.7270874847e+00,
 -1.7220652946e+00,
 -1.6602926834e+00,
 -1.6602926834e+00]
)



# these file names should still be good
dz_en_file = "HHep_fci_cc_pVDZ_Energies.npy"
dz_mu_file = "HHep_fci_cc_pVDZ_Dipoles.npy"
tz_en_file = "HHep_fci_cc_pVTZ_Energies.npy"
tz_mu_file = "HHep_fci_cc_pVTZ_Dipoles.npy"
qz_en_file = "HHep_fci_cc-pVQZ_Energies.npy"
qz_mu_file = "HHep_fci_cc-pVQZ_Dipoles.npy"


dz_en = np.load(dz_en_file)
dz_mu = np.load(dz_mu_file)
tz_en = np.load(tz_en_file)
tz_mu = np.load(tz_mu_file)
qz_en = np.load(qz_en_file)
qz_mu = np.load(qz_mu_file)

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
tz_inst = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)
qz_inst = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)

# number of adiabatic states for each basis set
N_el_dz = len(dz_en)
N_el_tz = len(tz_en)
N_el_qz = len(qz_en)
leng = 10

dz_dN = int(0.05 * N_el_dz)
tz_dN = int(0.05 * N_el_tz)
qz_dN = int(0.05 * N_el_qz)
print(F" Length of dz is {N_el_dz}, increment is {dz_dN}")
print(F" Length of tz is {N_el_tz}, increment is {tz_dN} total is {leng * tz_dN}")
print(F" Length of qz is {N_el_qz}, increment is {qz_dN}")


# array of number of electronic states in increments of 5% of total adiabatic states beginning from minimal basis
N_el_list_dz = np.linspace(dz_dN, leng * dz_dN, leng, dtype=int)
print(N_el_list_dz)
N_el_list_tz = np.linspace(tz_dN, leng * tz_dN, leng, dtype=int)
print(N_el_list_tz)
N_el_list_qz = np.linspace(qz_dN, leng * qz_dN, leng, dtype=int)
print(N_el_list_qz)
#N_el_list_tz = [2]
#N_el_list_dz = [2]


N_ph = 20
omega_dz = 0.9760568251
omega_tz = 0.9654959009
omega_qz = 0.9637811053



lambda_vector = np.array([0., 0., 0.02])

# double zeta arrays
dz_energies = np.zeros((leng, 10))
qed_fci_dz_energy = []

# triple zeta arrays
tz_energies = np.zeros((leng, 10))
qed_fci_tz_energy = []

# quadruple zeta arrays
qz_energies = np.zeros((leng, 10))
qed_fci_qz_energy = []

for i in range(leng):
    fast_start = time.time()
    # build double zeta Hamiltonian
    dz_inst.fast_build_pcqed_pf_hamiltonian(N_el_list_dz[i], N_ph, omega_dz, lambda_vector, dz_en, dz_mu, neglect_DSE=False)
    # build triple zeta Hamiltonian
    tz_inst.fast_build_pcqed_pf_hamiltonian(N_el_list_tz[i], N_ph, omega_tz, lambda_vector, tz_en, tz_mu, neglect_DSE=False)
    # build quadruple zeta Hamiltonian
    qz_inst.fast_build_pcqed_pf_hamiltonian(N_el_list_qz[i], N_ph, omega_qz, lambda_vector, qz_en, qz_mu, neglect_DSE=False)
    fast_end = time.time()
    dt = fast_end - fast_start
    print(F"Fast build took {dt} seconds")
    #print(dz_inst.PCQED_pf_eigs[0:10])
    # store double zeta energies
    dz_energies[i,:] = np.copy(dz_inst.PCQED_pf_eigs[:10])
    qed_fci_dz_energy.append(qed_fci_dz_large_lambda[0])
    # store triple zeta energies 
    tz_energies[i,:] = np.copy(tz_inst.PCQED_pf_eigs[:10])
    qed_fci_tz_energy.append(qed_fci_tz_large_lambda[0])
    # store quadruple zeta energies 
    qz_energies[i,:] = np.copy(qz_inst.PCQED_pf_eigs[:10])
    qed_fci_qz_energy.append(qed_fci_qz_large_lambda[0])
    
    #timing_list.append(dt)
