import numpy as np

def ground_disp_Mair_1993(x, z0, K, D, vl):
    """
    Calculate the 2D Gaussian based ground deformation. 
    The vertical displacement comes from Mair, R. J., Taylor, R. N., & Bracegirdle,
    A. (1993). Subsurface settlement profiles above tunnels in clays. Geotechnique,
    43(2), 315-320.
    The horizontal displacement comes from O'reilly, M. P., & New, B. M. (1982).
    Settlements above tunnels in the United Kingdom-their magnitude and prediction
    (No. Monograph).

    Parameters
    ----------
    x : float or np.array(dtype=float)
        Distance from tunnel center. Unit: m
    z0 : float or np.array(dtype=float)
        Tunnel center depth. Unit: m
    K : float or np.array(dtype=float)
        Trough width parameter. Unit: Unitless
    D : float
        Tunnel diameter. Unit: m
    vl : float or np.array(dtype=float)
        Tunnel volume loss. Unit: the actual quantity, not percentages

    Returns
    -------
    tuple
        A tuple containing:
        - The vertical displacement (float  or np.array(dtype=float)) in the unit
          of m and upward as positive value.
        - The horizontal displacement (float  or np.array(dtype=float)) in the unit
          of m and positive direction is the same with x.
    """
    # Eq. 3 of Mair
    i = K*z0
    # Eq. 6a of Mair
    s_max = 0.313 * vl * D**2 / i
    # EQ. 1 of Mair
    dv = - s_max * np.exp(-np.power(x, 2)/2/np.power(i,2))
    # Eq. 6 of O'Reilly
    dh = dv * x / z0
    # Return the vertical and horizontal disp
    return dv, dh


# def ground_disp_Mair_1993_regional(chainage, statis):
#     chainage