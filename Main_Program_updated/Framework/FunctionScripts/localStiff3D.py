import numpy as np
def moveResultLocation(model, d_a, numNodes):
    model.beam_DispL += model.beam_RotaT * d_a

    total_disp = np.zeros([numNodes * 6, 1])
    total_disp[0::6] = model.beam_DispL.reshape(-1, 1)  # [m] Longitudinal / axial (x)
    total_disp[1::6] = model.beam_DispT.reshape(-1, 1)  # [m] Transversal (y)
    total_disp[2::6] = model.beam_DispV.reshape(-1, 1)  # [m] Vertical (z)
    total_disp[3::6] = model.beam_RotaL.reshape(-1, 1)  # [rad] Horizontal rotations
    total_disp[4::6] = model.beam_RotaT.reshape(-1, 1)  # [rad] Transversal rotations
    total_disp[5::6] = model.beam_RotaV.reshape(-1, 1)  # [ras] Vertical rotations
    return model, total_disp


def localTim3Dold(EA, EIy, EIz, GAy, GAz, l, G, K):
    '''
    Calculates the local 3D stiffness matrix for a Timoshenko beam
    Args:
        EA: Axial stiffness
        EIy: Bending stiffness around the z axis (Moment rotates around z-axis)
        EIz: Bending stiffness around the y axis (Moment rotates around y-axis)
        GAy: Shear stiffness in the x-y plane
        GAz: Shear stiffness in the x-z plane
        l: Length of the beam element
        G: Isotropic Shear stiffness
        K: Saint Venant torsional stiffness

    Returns: K_local
    '''

    phi_y = 12 * (EIz / (GAy * l ** 2))
    phi_z = 12 * (EIy / (GAz * l ** 2))
    kz = EIz / ((1 + phi_y) * l ** 3)
    ky = EIy / ((1 + phi_z) * l ** 3)
    K_local = np.zeros([12, 12])

    # %% Indexing the axial stiffness
    K_local[0, 0] = EA / l
    K_local[0, 6] = -EA / l
    K_local[6, 0] = -EA / l
    K_local[6, 6] = EA / l

    # %% Indexing the top left quadrant
    K_local[1, 1] = 12 * kz  # kz11
    K_local[1, 5] = 6 * l * kz  # kz12

    K_local[2, 2] = 12 * ky  # ky11
    K_local[2, 4] = -6 * l * ky  # ky12

    K_local[3, 3] = G * K / l  # St. Venant torsion

    K_local[4, 2] = -6 * l * ky  # ky12
    K_local[4, 4] = (4 + phi_z) * (l ** 2) * ky  # ky22

    K_local[5, 1] = 6 * l * kz  # kz12
    K_local[5, 5] = (4 + phi_y) * (l ** 2) * kz  # kz12

    # %% Indexing the top right quadrant
    K_local[1, 7] = -12 * kz  # kz13
    K_local[1, 11] = 6 * l * kz  # kz14

    K_local[2, 8] = -12 * ky  # ky13
    K_local[2, 10] = -6 * l * ky  # ky14

    K_local[3, 9] = -G * K / l  # St. Venant torsion

    K_local[4, 8] = 6 * l * ky  # ky23
    K_local[4, 10] = (2 - phi_z) * (l ** 2) * ky  # ky24

    K_local[5, 7] = -6 * l * kz  # kz23
    K_local[5, 11] = (2 - phi_y) * (l ** 2) * kz  # kz24

    # %% Indexing the bottom left quadrant
    K_local[7, 1] = -12 * kz  # kz13
    K_local[7, 5] = -6 * l * kz  # kz23

    K_local[8, 2] = -12 * ky  # ky13
    K_local[8, 4] = 6 * l * ky  # ky23

    K_local[9, 3] = -G * K / l  # St. Venant torsion

    K_local[10, 2] = -6 * l * ky  # ky14
    K_local[10, 4] = (2 - phi_z) * (l ** 2) * ky  # ky24

    K_local[11, 1] = 6 * l * kz  # kz14
    K_local[11, 5] = (2 - phi_y) * (l ** 2) * kz  # kz24

    # %% Indexing the bottom right quadrant
    K_local[7, 7] = 12 * kz  # kz33
    K_local[7, 11] = -6 * l * kz  # kz34

    K_local[8, 8] = 12 * ky  # ky33
    K_local[8, 10] = 6 * l * ky  # ky34

    K_local[9, 9] = G * K / l  # St. Venant torsion

    K_local[10, 8] = 6 * l * ky  # ky34
    K_local[10, 10] = (4 + phi_z) * (l ** 2) * ky  # ky44

    K_local[11, 7] = -6 * l * kz  # kz34
    K_local[11, 11] = (4 + phi_y) * (l ** 2) * ky  # kz44

    return K_local


def compute_internal_forces(EA, EIy, EIz, GAy, GAz, l, G, K, total_disp, num_nodes):
    # Calculate the local stiffness matrix
    K_local = localTim3D(EA, EIy, EIz, GAy, GAz, l, G, K)
    # K_local = KBern3D_foot_TIM_dNA(1E6, 9, 0.5, l, 2.6, 0.3, 4.5)

    # Number of elements
    num_elements = num_nodes - 1

    # Initialize the internal forces vectors
    F_M_deltaT_el = np.zeros(2 * num_elements)
    F_N_deltaT_el = np.zeros(2 * num_elements)
    F_S_deltaT_el = np.zeros(2 * num_elements)

    # Compute the internal forces for each element
    for j in range(num_elements):
        # Extract the nodal displacements for the current element
        element_disp = np.zeros((12, 1))
        element_disp[0:6] = total_disp[j * 6:(j + 1) * 6]
        if (j + 1) * 6 < len(total_disp):
            element_disp[6:12] = total_disp[(j + 1) * 6:(j + 2) * 6]

        # Compute the internal forces for the current element
        Fin_P = np.dot(K_local, element_disp)

        # Store the internal forces in the vectors
        F_M_deltaT_el[2 * j] = -Fin_P[4, 0]
        F_M_deltaT_el[2 * j + 1] = Fin_P[10, 0]
        F_N_deltaT_el[2 * j] = -Fin_P[0, 0]
        F_N_deltaT_el[2 * j + 1] = Fin_P[6, 0]
        F_S_deltaT_el[2 * j] = -Fin_P[2, 0]
        F_S_deltaT_el[2 * j + 1] = Fin_P[8, 0]

    return F_M_deltaT_el, F_N_deltaT_el, F_S_deltaT_el


# Method used by Jinyan for comparison
def KBern3D_foot_TIM_dNA(E, d_foot, b_foot, dx, EGratio, ni_str, d_NA):
    '''
    Calculates the local 3D stiffness matrix for a Timoshenko beam with neutral axis distance d_NA
    Args:
        E: Young's modulus
        d_foot: Cross-sectional dimension
        b_foot: Cross-sectional dimension
        dx: Length of the beam element
        EGratio: Ratio of Young's modulus to shear modulus
        ni_str: Poisson's ratio
        d_NA: Neutral axis distance

    Returns: Kelem
    '''

    A = b_foot * d_foot
    I22 = b_foot * pow(d_foot, 3) / 12
    I33 = d_foot * pow(b_foot, 3) / 12
    I11 = I22 + I33
    L = dx
    k = 10 * (1 + ni_str) / (12 + 11 * ni_str)
    G = E / EGratio
    phi2 = 12 * E * I22 / k / G / A / (L * L)
    phi3 = 12 * E * I33 / k / G / A / (L * L)
    phi2Bar = 1 / (1 + phi2)
    phi3Bar = 1 / (1 + phi3)

    K = np.zeros([12, 12])
    K[0, 0] = E * A / L
    K[0, 4] = 1 * E * A / L * d_NA
    K[0, 6] = (-1) * E * A / L
    K[0, 10] = -1 * E * A / L * d_NA
    K[1, 1] = 12 * E * I33 * phi3Bar / pow(L, 3)
    K[1, 5] = 6 * E * I33 * phi3Bar / pow(L, 2)
    K[1, 7] = -12 * E * I33 * phi3Bar / pow(L, 3)
    K[1, 11] = 6 * E * I33 * phi3Bar / pow(L, 2)
    K[2, 2] = 12 * E * I22 / pow(L, 3) * phi2Bar
    K[2, 4] = -6 * E * I22 / pow(L, 2) * phi2Bar
    K[2, 8] = -12 * E * I22 / pow(L, 3) * phi2Bar
    K[2, 10] = -6 * E * I22 / pow(L, 2) * phi2Bar
    K[3, 3] = G * I11 / L
    K[3, 9] = (-1) * G * I11 / L
    K[4, 0] = E * A / L * d_NA
    K[4, 2] = -6 * E * I22 / pow(L, 2) * phi2Bar
    K[4, 4] = (4 + phi2) * E * I22 / L * phi2Bar + E * A / L * pow(d_NA, 2)
    K[4, 6] = - E * A / L * d_NA
    K[4, 8] = 6 * E * I22 / pow(L, 2) * phi2Bar
    K[4, 10] = (2 - phi2) * E * I22 / L * phi2Bar
    K[5, 1] = 6 * E * I33 * phi3Bar / pow(L, 2)
    K[5, 5] = (4 + phi3) * phi3Bar * E * I33 / L
    K[5, 7] = -6 * phi3Bar * E * I33 / pow(L, 2)
    K[5, 11] = (2 - phi3) * phi3Bar * E * I33 / L
    K[6, 0] = -1 * E * A / L
    K[6, 4] = - E * A / L * d_NA
    K[6, 6] = E * A / L
    K[6, 10] = E * A / L * d_NA
    K[7, 1] = -12 * E * I33 * phi3Bar / pow(L, 3)
    K[7, 5] = -6 * E * I33 * phi3Bar / pow(L, 2)
    K[7, 7] = 12 * E * I33 * phi3Bar / pow(L, 3)
    K[7, 11] = -6 * E * I33 * phi3Bar / pow(L, 2)
    K[8, 2] = -12 * E * I22 * phi2Bar / pow(L, 3)
    K[8, 4] = 6 * E * I22 / pow(L, 2) * phi2Bar
    K[8, 8] = 12 * E * I22 / pow(L, 3) * phi2Bar
    K[8, 10] = 6 * E * I22 / pow(L, 2) * phi2Bar
    K[9, 3] = -1 * G * I11 / L
    K[9, 9] = G * I11 / L
    K[10, 0] = - E * A / L * d_NA
    K[10, 2] = -6 * E * I22 / pow(L, 2) * phi2Bar
    K[10, 4] = (2 - phi2) * E * I22 / L * phi2Bar - E * A / L * pow(d_NA, 2)
    K[10, 6] = E * A / L * d_NA
    K[10, 8] = 6 * E * I22 / pow(L, 2) * phi2Bar
    K[10, 10] = (4 + phi2) * E * I22 / L * phi2Bar + E * A / L * pow(d_NA, 2)
    K[11, 1] = 6 * E * phi3Bar * I33 / pow(L, 2)
    K[11, 5] = (2 - phi3) * phi3Bar * E * I33 / L
    K[11, 7] = -6 * phi3Bar * E * I33 / pow(L, 2)
    K[11, 11] = (4 + phi3) * E * phi3Bar * I33 / L

    return K





def localTim3D(EA, EIy, EIz, GAy, GAz, l, G, J):
    '''
    Calculates the local 3D stiffness matrix for a Timoshenko beam
    Args:
        EA: Axial stiffness
        EIy: Bending stiffness around the z axis (Moment rotates around z-axis)
        EIz: Bending stiffness around the y axis (Moment rotates around y-axis)
        GAy: Shear stiffness in the x-y plane
        GAz: Shear stiffness in the x-z plane
        l: Length of the beam element
        G: Isotropic Shear stiffness
        K: Saint Venant torsional stiffness

    Returns: K_local
    '''

    # Due to convension, EIy = EIzz and likewise
    EIyy = EIz
    EIzz = EIy

    # New matrix:

    phi_z = 12 * EIzz / (GAz * (l ** 2))
    phi_z_bar = 1 / (1 + phi_z)

    phi_y = 12 * EIyy / (GAy * (l ** 2))
    phi_y_bar = 1 / (1 + phi_y)


    K_local = np.zeros([12, 12])

    # %% Indexing the axial stiffness
    K_local[0, 0] = EA / l
    K_local[0, 6] = -EA / l
    K_local[6, 0] = -EA / l
    K_local[6, 6] = EA / l


    # %% Indexing torsion
    K_local[3, 3] = G * J / l  # St. Venant torsion
    K_local[3, 9] = -G * J / l  # St. Venant torsion
    K_local[9, 3] = -G * J / l  # St. Venant torsion
    K_local[9, 9] = G * J / l  # St. Venant torsion


    # %% Indexing the 2nd row / col
    K_local[1, 1] = 12 * phi_z_bar * EIzz / l ** 3  # kz11
    K_local[1, 5] = 6 * phi_z_bar * EIzz / l ** 2  # kz12
    K_local[5, 1] = 6 * phi_z_bar * EIzz / l ** 2  # kz21

    K_local[7, 1] = -12 * phi_z_bar * EIzz / l ** 3
    K_local[1, 7] = -12 * phi_z_bar * EIzz / l ** 3

    K_local[11, 1] = 6 * phi_z_bar * EIzz / l ** 2
    K_local[1, 11] = 6 * phi_z_bar * EIzz / l ** 2

    # Now col 3
    K_local[2, 2] = 12 * phi_y_bar * EIyy / l ** 3  # ky11
    K_local[2, 4] = -6 * phi_y_bar * EIyy / l ** 2
    K_local[4, 2] = -6 * phi_y_bar * EIyy / l ** 2

    K_local[2, 8] = -12 * phi_y_bar * EIyy / l ** 3
    K_local[8, 2] = -12 * phi_y_bar * EIyy / l ** 3

    K_local[2, 10] = -6 * phi_y_bar * EIyy / l ** 2
    K_local[10, 2] = -6 * phi_y_bar * EIyy / l ** 2

    # Now Col 4 - Already done as this is torsion

    # Now Col 5
    K_local[4, 4] = (4 + phi_y) * phi_y_bar * EIyy / l  # ky11
    K_local[4, 8] = 6 * phi_y_bar * EIyy / l ** 2
    K_local[8, 4] = 6 * phi_y_bar * EIyy / l ** 2

    K_local[4, 10] = (2 - phi_y) * phi_y_bar * EIyy / l
    K_local[10, 4] = (2 - phi_y) * phi_y_bar * EIyy / l

    # Col 6
    K_local[5, 5] = (4 + phi_z) * phi_z_bar * EIzz / l  # ky11
    K_local[5, 7] = -6 * phi_z_bar * EIzz / l ** 2
    K_local[7, 5] = -6 * phi_z_bar * EIzz / l ** 2

    K_local[5, 11] = (2 - phi_z) * phi_z_bar * EIzz / l
    K_local[11, 5] = (2 - phi_z) * phi_z_bar * EIzz / l

    # Col 7 DOne (Axial)

    # Col 8
    K_local[7, 7] = 12 * phi_z_bar * EIzz / l ** 3
    K_local[7, 11] = -6 * phi_z_bar * EIzz / l ** 2
    K_local[11, 7] = -6 * phi_z_bar * EIzz / l ** 2

    # Col 9
    K_local[8, 8] = 12 * phi_y_bar * EIyy / l ** 3
    K_local[8, 10] = 6 * phi_y_bar * EIyy / l ** 2
    K_local[10, 8] = 6 * phi_y_bar * EIyy / l ** 2

    # Col 10 torsion skip

    # Col 11
    K_local[10, 10] = (4 + phi_y) * phi_y_bar * EIyy / l  # ky11

    # Col 12
    K_local[11, 11] = (4 + phi_z) * phi_z_bar * EIzz / l  # ky11


    '''
    K_local[2, 4] = -6 * phi_y_bar * EIyy / l ** 2  # ky12
    '''




    return K_local

