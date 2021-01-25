import time
import numpy as np
np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4

# Memory for Psi4 in GB
psi4.set_memory('2 GB')
psi4.core.set_output_file("test.dat", False)

# Memory for numpy in GB
numpy_memory = 2


def readXYZ(file):
    f = open(file, "r")
    lines = f.readlines()
    filelength = len(lines)
    progress = 0
    geomcount = 0
    geom = []
    while progress < filelength:
        tmpgeom = ""
        length = int(lines[progress])+2
        rangestart = progress
        rangeend = progress + length
        for i in range(rangestart, rangeend):
            tmpgeom = tmpgeom + lines[i]
        geom = geom + [tmpgeom]
        geomcount = geomcount + 1
        progress = progress + length
    f.close()
    return(''.join((geom)))

mol = psi4.geometry("""
O  0.0000   0.0000   -0.075791843589   
H  0.0000 -0.866811828967   0.601435779270   
H  0.00     0.866811828967   0.601435779270 
symmetry c1
no_reorient 
no_com
""")

psi4.set_options({'basis': '6-31g**',
                  'scf_type': 'direct',
                  'guess': 'core',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8,
                  'PERTURB_H': True,
                  'PERTURB_WITH':'DIPOLE',
                  'PERTURB_DIPOLE':[0.,0.,0.0001]})
scf_e, scf_wfn = psi4.energy('scf', return_wfn=True)


