import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator

# setup basic arguments to create an instance of the PFHamiltonianGenerator class
mol_str = """
H
He 1 0.776011
1 1
symmetry c1
"""


options_dict = {
    "basis": "cc-pVQZ",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
}


cavity_free_dict = {
    'omega_value' : 0,
    'lambda_vector' : np.array([0, 0, 0.0]),
    'ci_level' : 'fci',   
    'full_diagonalization' : True,
    'number_of_photons' : 0, 
}

# create the instance of our PFHamiltonianGenerator class
instance = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)
# figure out the number of states
_dim = len(instance.CIeigs)
# create array of states
states = np.linspace(0, _dim-1, _dim, dtype=int)
# compute the dipole moment for all states
mu_array = instance.compute_dipole_moments(states)

file_string = "HHep_fci_cc_pVQZ"
E_string = file_string + "_Energies"
Mu_string = file_string + "_Dipoles"

#np.save(C_string, test_pf.CIvecs)
np.save(E_string, instance.CIeigs)
np.save(Mu_string, mu_array)                      
