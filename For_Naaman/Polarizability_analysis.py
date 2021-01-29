import numpy as np
from numpy import linalg as LA

HS_25_polarizability = np.array([[ -616.06509371,  -338.23860565,   168.61275949],
 [ -339.53209953, -3271.1258288,     49.97796199],
 [  169.52823124,    49.92681813, -3248.46416257]])
HS_13_polarizability = np.array( [[ -632.12055991,   361.8850252,   -105.08752708],
 [  362.5104567,  -2561.1274773,    -86.79063915],
 [ -105.24760956,   -86.84916529, -2747.83281691]])
HS_4_polarizability = np.array( [[ -835.28272287,     9.73158473,     5.52674842],
 [    9.73998582, -1875.59484557,     2.56536624],
 [    5.51324309,     2.56795261, -1894.69617352]])
HS_40_polarizability = np.array( [[ -862.9462555,    -51.67449766,    13.62415435],
 [  -51.70552815, -2339.69793078,     6.70147548],
 [   13.61809558,     6.69942881, -2349.18932917]])


def normalize_rows(x: np.ndarray):
    """
    function that normalizes each row of the matrix x to have unit length.

    Args:
     ``x``: A numpy matrix of shape (n, m)

    Returns:
     ``x``: The normalized (by row) numpy matrix.
    """
    return x/LA.norm(x, ord=2, axis=1, keepdims=True)


def polarizability_analisys_norm(Pol):
    _Pol = normalize_rows(Pol)
    xx, yy, zz = _Pol.diagonal()
    isotropy = (xx+yy+zz)/3
    anisotropy = (0.5 * ((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2))**0.5
    print('Isotropic polarizability: {:}'.format(isotropy))
    print('Polarizability anisotropy: {:}'.format(anisotropy))
    print('Nromalized dipole polarizability tensor: {:}'.format(_Pol))


if __name__ == '__main__':
    print("RESULTS IN DECREASING ANISOTROPY")
    print("================================")
    print("HS_25")
    print("----")
    polarizability_analisys_norm(HS_25_polarizability)
    print("HS_13")
    print("-----")
    polarizability_analisys_norm(HS_13_polarizability)
    print("HS_40")
    print("-----")
    polarizability_analisys_norm(HS_40_polarizability)
    print("HS_4")
    print("----")
    polarizability_analisys_norm(HS_4_polarizability)
