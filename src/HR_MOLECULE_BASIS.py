import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator
import scipy

#Molecule Geometry
mol_str = """
H
H 1 0.74
symmetry c1
"""

# options for the PFHamiltonian Generator class - include cavity effects
cavity_dict = {
    'omega_value' : 0.7320720628787358, #Change to guess for eq bond length
    'lambda_vector' : np.array([0, 0, 0.000]),
    'ci_level' : 'fci',
    'full_diagonalization' : True,
    'number_of_photons' : 1, #<== this is a minimal photon basis, should explore increasing this
}

# options for PFHamiltonian Generator class - exclude cavity effects
cavity_free_dict = {
    'omega_value' : 0.0,
    'lambda_vector' : np.array([0, 0, 0.0]),
    'ci_level' : 'fci',
    'full_diagonalization' : True,
    'number_of_photons' : 0,
}

#Change to atoms of interest
mol_tmpl = """
H
H 1 **R**
symmetry c1
"""
options_dict = {
    "basis": "6-311G",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
    'num_roots' : 2
}

# number of bondlengths in the scan
N_R = 100

# number of electronic states to save
N_el = 8

# array  for energies inside the cavity
cavity_E_array = np.zeros((N_R, N_el))

r_data = np.linspace(0.5, 2.0, N_R)
psi4.set_options(options_dict)
fci_S0 = []
fci_S1 = []
r_idx = 0
for r in r_data:
    mol_str = mol_tmpl.replace("**R**", str(r))
    mol = psi4.geometry(mol_str)
    scf_e, wfn = psi4.energy('SCF', return_wfn=True)
    fci_energy, wfn = psi4.energy('fci',ref_wfn=wfn, return_wfn=True)
    fci_S0.append(wfn.variable("CI ROOT 0 TOTAL ENERGY"))
    fci_S1.append(wfn.variable("CI ROOT 1 TOTAL ENERGY"))
    cav = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
    cavity_E_array[r_idx,:] = cav.CIeigs[:N_el]
    r_idx += 1

fci_S0_array = np.array(fci_S0)
fci_S1_array = np.array(fci_S1)
#For every line, change MOLECULE to the molecule of interest (and basis set if not at 6-311G)
np.save("fci_S0_MOLECULE_6311G", fci_S0) 
np.save("fci_S1_MOLECULE_6311G", fci_S1)
np.save("fci_cavity_arrays_MOLECULE_6311G", cavity_E_array)
np.save("fci_r_array_MOLECULE_6311G", r_data)

print('Do it fart?')





