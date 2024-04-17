import numpy as np
import sys
import psi4
from helper_PFCI import PFHamiltonianGenerator
import scipy

num_roots = 12

#Molecule Geometry
mol_str = """
Li
H 1 1.068
symmetry c1
"""

#/Users/proden/Code/SCQED-PCQED/src/HR_LIH_6311g.py


# options for PFHamiltonian Generator class - exclude cavity effects
cavity_dict = {
    'omega_value' : 0.105,
    'lambda_vector' : np.array([0, 0, 0.050]),
    'ci_level' : 'fci',
    'number_of_photons' : 5,
    'davidson_roots' : num_roots,
    'photon_number_basis' : False,
    'canonical_mos':True,
    'coherent_state_basis': True
}

# options for PFHamiltonian Generator class - exclude cavity effects
cavity_free_dict = {
    'omega_value' : 0.0,
    'lambda_vector' : np.array([0, 0, 0.0]),
    'ci_level' : 'fci',
    'number_of_photons' : 0,
    'davidson_roots' : num_roots,
    'photon_number_basis' : True,
    'canonical_mos':True,
    'coherent_state_basis': False
}

#Change to atoms of interest
mol_tmpl = """
Li
H 1 **R**
symmetry c1
"""

# number of bondlengths in the scan
N_R = 1000


options_dict = {
    "basis": "6-311g",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
    'num_roots' : num_roots
}


#determine how many states to sve for cavity free
cav_free = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)
N_el_cav_free = cav_free.CIeigs.shape[0]

cav= PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
N_el_cav = cav.CIeigs.shape[0]



options_dict = {
    "basis": "6-311g",
    "scf_type": "pk",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
    'num_roots' : num_roots
}


# array  for energies inside the cavity
cavity_E_array = np.zeros((N_R, N_el_cav))
cavity_free_E_array = np.zeros((N_R, N_el_cav_free))

r_data = np.linspace(1.0, 3, N_R)
psi4.set_options(options_dict)


r_idx = 0




for r in r_data:
    mol_str = mol_tmpl.replace("**R**", str(r))
    cav_free = PFHamiltonianGenerator(mol_str, options_dict, cavity_free_dict)
    cavity_free_E_array[r_idx,:] = cav_free.CIeigs[:N_el_cav_free]
    r_idx += 1

 #For every line, change MOLECULE to the molecule of interest (and basis set if not at 6-311G)
np.save("cs_5_ph_om_0_105_fci_cavity_free_array_LIH_6311g", cavity_free_E_array)

r_idx = 0







lambda_values = [0.005, 0.01, 0.015,0.020, 0.025,0.03,0.04,0.05,0.06,0.08,0.1]
for i in range(0,len(lambda_values)):
 r_idx = 0
 cavity_E_array = np.zeros((N_R, N_el_cav))
 Dipoles = np.zeros((num_roots, num_roots, 3, N_R))
 cavity_dict = {
    'omega_value' : 0.105,
    'lambda_vector' : np.array([0, 0, lambda_values[i]]),
    'ci_level' : 'fci',
    'number_of_photons' : 5,
    'davidson_roots' : num_roots,
    'photon_number_basis' : False,
    'canonical_mos':True,
    'coherent_state_basis':True
 }

 for r in r_data:
    mol_str = mol_tmpl.replace("**R**", str(r))

    cav = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
    cavity_E_array[r_idx,:] = cav.CIeigs[:N_el_cav]
    Dipoles[:,:,:,r_idx] = np.copy(cav.dipole_array)
    r_idx += 1

 #For every line, change MOLECULE to the molecule of interest (and basis set if not at 6-311G) 
 np.save("cs_5_ph_om_0_105_fci_cavity_arrays_LIH_6311g"+ str(lambda_values[i]).replace(".","_"), cavity_E_array)
 np.save("cs_5_ph_om_0_105_fci_dipoles_LIH_6311g"+str(lambda_values[i]).replace(".","_"), Dipoles)
np.save("cs_5_ph_om_0_105_fci_r_array_LIH_6311G", r_data)

