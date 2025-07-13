import sys
sys.path.append('FunctionScripts')
from localStiff3D import *



# Validation of inputs into stiffness matrix:
E = 1000
nu = 0.3
G = E / (2 + 2 * nu)


L = 0.2
h = 1
b = 1



Izz = (1 / 12) * (h ** 3 * b) * 2
Iyy = (1 / 12) * (b ** 3 * h)
A = h * b
As = (10 + 10 * nu) / (11 * nu + 12) * A


phi_z = 12 * E * Izz / (G * As * (L ** 2))
phi_z_bar = 1 / (1 + phi_z)

phi_y = 12 * E * Iyy / (G * As * (L ** 2))
phi_y_bar = 1 / (1 + phi_y)




K_old = localTim3D(E*A, E*Izz, E*Iyy, G*As, G*As, L, G, 0)
# print('OldStiffness matrix', K_old)
# '''
# From paper: K[0,0] = EA / L
print('stiffness matrix result = ', K_old[0, 0], ', result form paper = ', E*A / L)

# From paper: K[1,1] = 12 * phi_z_bar * E*Izz / L^3
print('stiffness matrix result = ', K_old[1, 1], ', result form paper = ', 12 * phi_z_bar * E * Izz / L ** 3)

# From paper: K[4,4] = ( 4 + phi_y) * phi_y_bar * EIyy / L
print('stiffness matrix result = ', K_old[4, 4], ', result form paper = ', (4 + phi_y) * phi_y_bar * E * Iyy / L)


# From paper: K[5,5] = ( 4 + phi_z) * phi_z_bar * EIzz / L
print('stiffness matrix result = ', K_old[5, 5], ', result form paper = ', (4 + phi_z) * phi_z_bar * E * Izz / L)


# From paper: K[8,4] = 6 * phi_y_bar * E * Iyy / L ** 2
print('stiffness matrix result = ', K_old[8, 4], ', result form paper = ', 6 * phi_y_bar * E * Iyy / L ** 2)
# '''


import numpy as np


def is_symmetric(matrix, tol=1e-8):
    return np.allclose(matrix, matrix.T, atol=tol)

print('Is the matrix symmetrical? Answer = ', is_symmetric(K_old))
