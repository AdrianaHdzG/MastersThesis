import numpy as np

def prepare_greenfield(z0, vlt, Kt, d, l_b, eccentricity, num):
    '''
    Computes greenfield vertical and horizontal displacements along the building
    due to tunneling-induced ground movements using a Gaussian settlement trough.

    Parameters:
    -----------
    z0 : float
        Depth from ground surface to tunnel axis [m].
    vlt : float
        Volume loss due to tunneling [%].
    Kt : float
        Trough width parameter (empirical constant).
    d : float
        Tunnel diameter [m].
    l_b : float
        Length of the building [m].
    eccentricity : float
        Horizontal offset of the building center from the tunnel centerline [-].
    num : int
        Number of points/nodes along the building.

    Returns:
    --------
    Uz : np.ndarray
        Vertical displacement profile along the building [m].
    Ux : np.ndarray
        Horizontal displacement profile along the building [m].
    x_coord : np.ndarray
        X-coordinates of the building domain [m].

    Notes:
    ------
    - The vertical displacement is computed using a Gaussian distribution.
    - The horizontal displacement is derived from the vertical profile using
      the assumption of volume conservation and geometry of the settlement trough.
    '''
    # Describe the domain
    domainStart = -(0.5 + eccentricity) * l_b
    domainEnd = (0.5 - eccentricity) * l_b

    x_coord = np.linspace(domainStart, domainEnd, num)  # Create domain vector

    i = Kt * z0  # Calculate horizontal distance to inflection point from tunnel centre line

    Uz_max = 1.25 * (d / 2) ** 2 * vlt / i  # Calculate maximum settlement

    Uz = np.zeros([num])
    Ux = np.zeros([num])
    for j in range(num):
        Uz[j] = -Uz_max * np.exp(-(x_coord[j] ** 2) / (2 * i ** 2))  # Vertical settlement
        Ux[j] = Uz[j] * x_coord[j] / z0  # Horizontal displacement

    return Uz, Ux, x_coord



