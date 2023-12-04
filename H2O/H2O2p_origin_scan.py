import psi4
from helper_PFCI import PFHamiltonianGenerator
import numpy as np

file_string = "H2O2p_origin_scan_631g_fci"


mol_tmpl = """
    O            0.000000000000     0.000000000000    **ZO**   
    H            0.000000000000    -0.790689573744    **ZH**
    H            0.000000000000     0.790689573744    **ZH**     
    2 1
    symmetry c1
    no_reorient
    no_com
"""


# origin z-coordinate of O and H atoms
O_origin =  -0.068516219320
H_origin = 0.543701060715


# number of displacements in scan
N_R = 21

# evenly spaced grid of displacements
d_array = np.linspace(0, 20, N_R)


# number of CI roots to compute 
num_roots = 20

# size of energy array will be num_roots x N_R
energy_array = np.zeros((num_roots, N_R))

# size of dipole array will be (num_roots x num_roots x 3 x N_R)
dipole_array = np.zeros((num_roots, num_roots, 3, N_R))


options_dict = {
        "basis": "6-31g",
        "scf_type": "pk",
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
}

cavity_dict = {
        'omega_value' : 0,
        'lambda_vector' : np.array([0, 0, 0.0]),
        'ci_level' : 'fci',
        'davidson_roots' : num_roots,
        'number_of_photons' : 0,
        'photon_number_basis' : True,
        'canonical_mos' : True,
        'coherent_state_basis' : False
}


# loop over each geometry, run fci calculation, add eigenvalues and dipole matrix elements
# to the appropriate slice of the energy_array and dipole_array
ctr = 0
for d in d_array:
    O_z_coord = O_origin + d
    H_z_coord = H_origin + d
    tmp = mol_tmpl.replace("**ZO**", str(O_z_coord))
    mol_str = tmp.replace("**ZH**", str(H_z_coord))
    mol = psi4.geometry(mol_str)
    print(mol_str)
    print(cavity_dict)
    print(options_dict)

    test_pf = PFHamiltonianGenerator(mol_str, options_dict, cavity_dict)
    energy_array[:,ctr] = np.copy(test_pf.CIeigs)
    dipole_array[:,:,:,ctr] = np.copy(test_pf.dipole_array)
    ctr += 1



# save eigenvalues to 
E_string = file_string + "_Energies"
Mu_string = file_string + "_Dipoles"
np.save(E_string, energy_array)
np.save(Mu_string, dipole_array)

