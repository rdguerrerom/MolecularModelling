import psi4
molecule = psi4.geometry(
    """symmetry c1 # this is important; some perturbations will break symmetry
    no_reorient
    O          -0.947809457408    -0.132934425181     0.000000000000
    H          -1.513924046286     1.610489987673     0.000000000000
    F           0.878279174340     0.026485523618     0.000000000000""")

pert = 0.001
lambdas = [pert, -pert, 2.0*pert, -2.0*pert]

#psi4.set_options({'basis': 'cc-pVDZ',
#                  'scf_type': 'direct',
#                  'guess': 'core',
#                  'e_convergence': 1e-8,
#                  'd_convergence': 1e-8})
psi4.set_options({'basis': '6-31g**',
                  'scf_type': 'direct',
                  'guess': 'core',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8,
                  'PERTURB_H': True,
                  'PERTURB_WITH':'DIPOLE'})
scf_e, scf_wfn = psi4.energy('MP2', return_wfn=True)

perturbed_energies = np.zeros((len(lambdas),3))
x_perturbed_dipoles = np.zeros((len(lambdas),3))
y_perturbed_dipoles = np.zeros((len(lambdas),3))
z_perturbed_dipoles = np.zeros((len(lambdas),3))

# start with a reference dipole calculation
properties(method, properties=['dipole'])
mu_x = psi4.core.variable(method + ' DIPOLE X') / psi_dipmom_au2debye
mu_y = psi4.core.variable(method + ' DIPOLE Y') / psi_dipmom_au2debye
mu_z = psi4.core.variable(method + ' DIPOLE Z') / psi_dipmom_au2debye
analytic_dipole = np.array([mu_x, mu_y, mu_z])

# now compute with different applied fields
for step, l in enumerate(lambdas):
    #set perturb_h true
    #set perturb_with dipole
    # x pertubation
    #set perturb_dipole [$l, 0, 0]
    psi4.set_options({'PERTURB_DIPOLE':[0.,0.,0.0001]})
    perturbed_energies[step,0] = properties(method, properties=['dipole'])
    x_perturbed_dipoles[step,0] = psi4.core.variable(method + ' DIPOLE X')
    x_perturbed_dipoles[step,1] = psi4.core.variable(method + ' DIPOLE Y')
    x_perturbed_dipoles[step,2] = psi4.core.variable(method + ' DIPOLE Z')
    # y pertubation
    set perturb_dipole [0, $l, 0]
    perturbed_energies[step,1] = properties(method, properties=['dipole'])
    y_perturbed_dipoles[step,0] = psi4.core.variable(method + ' DIPOLE X')
    y_perturbed_dipoles[step,1] = psi4.core.variable(method + ' DIPOLE Y')
    y_perturbed_dipoles[step,2] = psi4.core.variable(method + ' DIPOLE Z')
    # z pertubation
    set perturb_dipole [0, 0, $l]
    perturbed_energies[step,2] = properties(method, properties=['dipole'])
    z_perturbed_dipoles[step,0] = psi4.core.variable(method + ' DIPOLE X')
    z_perturbed_dipoles[step,1] = psi4.core.variable(method + ' DIPOLE Y')
    z_perturbed_dipoles[step,2] = psi4.core.variable(method + ' DIPOLE Z')

# use 3- and 5-point finite difference formulae to compute the dipole
mu_3pt = (perturbed_energies[0] - perturbed_energies[1]) / (2.0*pert)
mu_5pt = (8.0*perturbed_energies[0] - 8.0*perturbed_energies[1] - perturbed_energies[2] + perturbed_energies[3]) / (12.0*pert)
# and the polarizability
a_x_3pt = (x_perturbed_dipoles[0] - x_perturbed_dipoles[1]) / (2.0*pert)
a_y_3pt = (y_perturbed_dipoles[0] - y_perturbed_dipoles[1]) / (2.0*pert)
a_z_3pt = (z_perturbed_dipoles[0] - z_perturbed_dipoles[1]) / (2.0*pert)
alpha_3pt = np.vstack((a_x_3pt, a_y_3pt, a_z_3pt))
a_x_5pt = (8.0*x_perturbed_dipoles[0] - 8.0*x_perturbed_dipoles[1] - x_perturbed_dipoles[2] + x_perturbed_dipoles[3]) / (12.0*pert)
a_y_5pt = (8.0*y_perturbed_dipoles[0] - 8.0*y_perturbed_dipoles[1] - y_perturbed_dipoles[2] + y_perturbed_dipoles[3]) / (12.0*pert)
a_z_5pt = (8.0*z_perturbed_dipoles[0] - 8.0*z_perturbed_dipoles[1] - z_perturbed_dipoles[2] + z_perturbed_dipoles[3]) / (12.0*pert)
alpha_5pt = np.vstack((a_x_5pt, a_y_5pt, a_z_5pt))

# let's see what happened!
np.set_printoptions(suppress=True)
print('Analytic dipole\n', analytic_dipole)
print('3point finite difference dipole\n', mu_3pt)
print('5point finite difference dipole\n', mu_5pt)

print('3point finite difference polarizability\n', alpha_3pt)
print('5point finite difference polarizability\n', alpha_5pt)
