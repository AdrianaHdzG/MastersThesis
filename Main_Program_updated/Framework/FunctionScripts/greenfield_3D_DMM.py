# FunctionScripts/greenfield_3D_DMM.py
import numpy as np
from math import erf, erfc

# ---- helpers (direct ports, minimal edits) ----
def _depth_parabolic(z_vec, Hw):
    z = np.asarray(z_vec).reshape(-1, 1)
    delta_w = (-6.0 / Hw) * z**2 + 6.0 * z
    return delta_w

def _depth_dmm(z_vec, Hw, C1, C2, C3):
    z = np.asarray(z_vec).reshape(-1, 1)
    delta_w = (
        C1 * (2*Hw - 2*z)
        + C2 * ((-6.0 / Hw) * z**2 + 6.0 * z)
        + C3 * (2.0 * z)
    )
    return delta_w

def _long_none(s_vec, Lwall, Hw, He_Hwratio):
    return np.ones((1, len(s_vec)), dtype=float)

def _long_mu(s_vec, Lwall, Hw, He_Hwratio):
    s = np.asarray(s_vec).reshape(1, -1)
    # width from Mu & Huang (kept as in your MATLAB)
    W = (Lwall/2.0) * (0.069*np.log((Hw/He_Hwratio)/Lwall) + 1.03)
    W = max(W, 0.05*Lwall)
    return np.exp(-np.pi * (s / W) ** 2)

def _long_roboski(s_vec, Lwall, Hw, He_Hwratio):
    s = np.asarray(s_vec).reshape(1, -1)
    a = 2.8
    c = 0.015 + 0.035*np.log(He_Hwratio*Hw/Lwall)
    d = 0.15  + 0.035*np.log(He_Hwratio*Hw/Lwall)
    num = a * ((Lwall/2.0 - np.abs(s)) + Lwall*c)
    den = max((0.5*Lwall - Lwall*d), 1e-9)
    return 1.0 - 0.5 * (1.0 - erfc(num/den))  # 1 - 0.5*erfc(x) = 0.5*(1+erf(x))

# Pinto/AF shaft solution: displacement from one equivalent spherical cavity
def _eq_shaft_3d_af(x, y, z, h, nu, a):
    R1 = np.sqrt(x**2 + y**2 + (z - h)**2)
    R2 = np.sqrt(x**2 + y**2 + (z + h)**2)

    fx = x * (1.0/R1**3 + (3 - 4*nu)/R2**3 - 6*z*(z + h)/R2**5)
    fy = y * (1.0/R1**3 + (3 - 4*nu)/R2**3 - 6*z*(z + h)/R2**5)
    fz = -1.0 * ((z - h)/R1**3 + 2*z/R2**3 - (3 - 4*nu)*(z + h)/R2**3 - 6*z*(z + h)**2/R2**5)

    # u = (1/3) * a^3 * f ; sign convention as in MATLAB (ux,uy negated)
    coeff = (a**3) / 3.0
    u_x = -coeff * fx
    u_y = -coeff * fy
    u_z =  coeff * fz
    return u_x, u_y, u_z

def run_greenfield_3D_line(
    *,
    # geometry
    Hw=19.0, L_x=9.5, L_y=32.0, He_Hwratio=1.0,
    # soil
    nu=0.499,
    # wall shape
    switch_shape=5, C1=0.0, C2=1.0, C3=None,
    # betas per wall
    beta_CCS_wall_1=0.075/100, beta_CCS_wall_2=0.075/100,
    beta_CCS_wall_3=0.075/100, beta_CCS_wall_4=0.075/100,
    # discretizations
    delta_z_cavities=None, delta_xyperimeter_cavities=2.5,
    # solution (1=single wall mirrored, 2=four walls analytical, 3=four walls semi-analytical)
    switch_solution_type=3,
    # output line (your building coordinates)
    building_offset=11.0, length_beam=12.0, num_nodes_3D=101, y0=0.0, z0=0.0,
):
    """
    Compute 3D greenfield displacements along a 1D building line (x-direction),
    at fixed y=y0 and z=z0. Returns the 2D-style outputs:
    horizontal_displacement, vertical_displacement, cavity_depth_vec, delta_wall_vector
    """
    if C3 is None:
        C3 = 1.0 - C1 - C2
    if abs((C1+C2+C3) - 1.0) > 1e-6 and switch_shape in (5, 50, 51):
        raise ValueError("C1+C2+C3 must be 1 for switch_shape=5/50/51")

    if delta_z_cavities is None:
        delta_z_cavities = 1.0  # default 1m vertical bins

    # ---- build the evaluation line (global coords) ----
    x_line = (L_x + building_offset) + np.linspace(0.0, length_beam, num_nodes_3D)
    y_line = np.full_like(x_line, y0, dtype=float)
    z_line = np.full_like(x_line, z0, dtype=float)

    # --- Taper weights along the evaluation line (SPy) ---
    # Mirrors your MATLAB masks for solution_type==3; for others, weights=1
    if switch_solution_type == 3:
        # Initialize as ones
        F1_line = np.ones_like(x_line)
        F2_line = np.ones_like(x_line)
        F3_line = np.ones_like(x_line)
        F4_line = np.ones_like(x_line)

        # Wall 1 (bottom, Y=-L_y): 1 at y=-Ly → 0 at y=+Ly (constant along x)
        if y0 >= -L_y:
            F1_line = np.maximum(0.0, 0.5 - 0.5*(y0 / L_y)) * np.ones_like(x_line)

        # Wall 2 (right, X=+L_x): 1 at x=+Lx → 0 at x=-Lx (varies with x for x<=+Lx; else stays 1)
        F2_line = np.where(
            x_line <= +L_x,
            np.maximum(0.0, 0.5 + 0.5*(x_line / L_x)),
            1.0
        )

        # Wall 3 (top, Y=+L_y): 1 at y=+Ly → 0 at y=-Ly (constant along x)
        if y0 <= +L_y:
            F3_line = np.maximum(0.0, 0.5 + 0.5*(y0 / L_y)) * np.ones_like(x_line)

        # Wall 4 (left, X=-L_x): 1 at x=-Lx → 0 at x=+Lx (varies with x for x>=-Lx; else stays 1)
        F4_line = np.where(
            x_line >= -L_x,
            np.maximum(0.0, 0.5 - 0.5*(x_line / L_x)),
            1.0
        )
    else:
        # For solution types 1 and 2, no taper (all ones)
        F1_line = F2_line = F3_line = F4_line = np.ones_like(x_line)
    
        
    # ---- wall perimeter midpoints (four sides of the box) ----
    # choose segment size; compute integer counts that exactly tile the lengths
    n1 = max(1, int(round(2*L_x / delta_xyperimeter_cavities)))
    n2 = max(1, int(round(2*L_y / delta_xyperimeter_cavities)))
    dx1 = (2*L_x) / n1
    dx2 = (2*L_y) / n2

    # wall 1 (bottom, y=-L_y): x mids
    xmid_w1 = np.linspace(-L_x + dx1/2, L_x - dx1/2, n1); ymid_w1 = -L_y*np.ones(n1)
    # wall 2 (right, x=+L_x): y mids
    ymid_w2 = np.linspace(-L_y + dx2/2, L_y - dx2/2, n2); xmid_w2 =  L_x*np.ones(n2)
    # wall 3 (top, y=+L_y): x mids (reversed traversal ok)
    xmid_w3 = np.linspace( L_x - dx1/2, -L_x + dx1/2, n1); ymid_w3 =  L_y*np.ones(n1)
    # wall 4 (left, x=-L_x): y mids
    ymid_w4 = np.linspace( L_y - dx2/2, -L_y + dx2/2, n2); xmid_w4 = -L_x*np.ones(n2)

    Xc_vec = np.concatenate([xmid_w1, xmid_w2, xmid_w3, xmid_w4])
    Yc_vec = np.concatenate([ymid_w1, ymid_w2, ymid_w3, ymid_w4])

    # index ranges per wall
    idx_w1 = np.arange(0, n1)
    idx_w2 = np.arange(n1, n1+n2)
    idx_w3 = np.arange(n1+n2, 2*n1+n2)
    idx_w4 = np.arange(2*n1+n2, 2*n1+2*n2)

    Lwall_1 = 2*L_x;  Lwall_2 = 2*L_y;  Lwall_3 = 2*L_x;  Lwall_4 = 2*L_y

    # ---- z discretization: cavities at mid-intervals ----
    z_wall_discr = np.arange(0.0, Hw + delta_z_cavities, delta_z_cavities)            # nodes
    z_cav        = (z_wall_discr[:-1] + z_wall_discr[1:]) * 0.5                       # mids

    # ---- choose depth and longitudinal functions ----
    if switch_shape in (3, 30, 31):
        depth_fun = lambda z: _depth_parabolic(z, Hw)
    elif switch_shape in (5, 50, 51):
        depth_fun = lambda z: _depth_dmm(z, Hw, C1, C2, C3)
    else:
        raise ValueError(f"Unknown switch_shape={switch_shape}")

    if switch_shape in (3, 5):
        long_fun = _long_none
    elif switch_shape in (30, 50):
        long_fun = _long_mu
    elif switch_shape in (31, 51):
        long_fun = _long_roboski
    else:
        raise ValueError(f"Unknown switch_shape={switch_shape}")

    # depth vectors
    delta_wall = depth_fun(z_wall_discr.reshape(-1, 1))   # Nd_w × 1
    delta_cav  = depth_fun(z_cav.reshape(-1, 1))          # Nd_c × 1

    # along-wall coordinate (centered by construction)
    s_w1 = xmid_w1; s_w2 = ymid_w2; s_w3 = xmid_w3; s_w4 = ymid_w4
    R1 = long_fun(s_w1, Lwall_1, Hw, He_Hwratio)  # 1 × Ns1
    R2 = long_fun(s_w2, Lwall_2, Hw, He_Hwratio)
    R3 = long_fun(s_w3, Lwall_3, Hw, He_Hwratio)
    R4 = long_fun(s_w4, Lwall_4, Hw, He_Hwratio)

    # outer-products: depth × along-wall
    DELTA_cav_1 = beta_CCS_wall_1 * (delta_cav @ R1)
    DELTA_cav_2 = beta_CCS_wall_2 * (delta_cav @ R2)
    DELTA_cav_3 = beta_CCS_wall_3 * (delta_cav @ R3)
    DELTA_cav_4 = beta_CCS_wall_4 * (delta_cav @ R4)
    DELTA_cav_all = np.concatenate([DELTA_cav_1, DELTA_cav_2, DELTA_cav_3, DELTA_cav_4], axis=1)

    # ---- symmetry factor ----
    if switch_solution_type == 1:
        m_sym = 2.0
    elif switch_solution_type == 2:
        m_sym = 1.0
    elif switch_solution_type == 3:
        m_sym = 2.0
    else:
        raise ValueError("switch_solution_type must be 1,2,or 3")

    # segment lengths along walls
    seg_len = np.concatenate([
        np.full(n1, dx1),
        np.full(n2, dx2),
        np.full(n1, dx1),
        np.full(n2, dx2),
    ])[None, :]  # shape (1, Ns)

    # equivalent cavity radii from volume elements
    Nd, Ns = DELTA_cav_all.shape
    vol_elem = m_sym * DELTA_cav_all * (delta_z_cavities * seg_len)  # Nd×Ns
    vol_elem = np.maximum(vol_elem, 0.0)
    a_cavity = ((3.0 / (4.0*np.pi)) * vol_elem) ** (1.0/3.0)         # Nd×Ns

    # ---- accumulate along the building line (y=y0,z=z0) ----
    Sx = np.zeros_like(x_line, dtype=float)
    Sy = np.zeros_like(x_line, dtype=float)  # unused
    Sz = np.zeros_like(x_line, dtype=float)

    # we only need ux, uz on the line → set y=0 by construction using translation
    # (the line's y0 is already baked in y_line, we subtract each stack y)
    index_sets = [idx_w1, idx_w2, idx_w3, idx_w4]
    fw_map = {1: F1_line, 2: F2_line, 3: F3_line, 4: F4_line}

    for w, idxs in enumerate(index_sets, start=1):
        Fw_line = fw_map[w]
        for j in idxs:
            dxj = Xc_vec[j]
            dyj = Yc_vec[j]

            # line translated to this stack
            xloc = x_line - dxj
            yloc = y_line - dyj
            zloc = z_line

            # sum depth contributions
            for i in range(Nd):
                a = a_cavity[i, j]
                if a <= 0.0:
                    continue
                h = z_cav[i]
                ux, uy, uz = _eq_shaft_3d_af(xloc, yloc, zloc, h, nu, a)
                Sx += Fw_line * ux
                Sy += Fw_line * uy
                Sz += Fw_line * uz

    # ---- mask inside the station (greenfield outside only) ----
    inside = (x_line >= -L_x) & (x_line <= L_x) & (np.abs(y_line) <= L_y)
    Sx = Sx * (~inside)
    Sy = Sy * (~inside)
    Sz = Sz * (~inside)

    # ---- produce 2D-style outputs ----
    horizontal_displacement = Sx.reshape(1, -1)   # (1, num_nodes) to match 2D shape
    transverse_displacement = Sy.reshape(1, -1)
    vertical_displacement   = Sz.reshape(1, -1)

    # For “delta_wall_vector” we return the vertical profile used for **one wall**
    # (parabolic/DMM without longitudinal degradation). This matches your 2D “dw vs depth” role.
    depth_vec = np.arange(0.0, Hw + (Hw/100.0), (Hw/100.0))
    delta_wall_vector = _depth_dmm(depth_vec, Hw, C1, C2, 1.0 - C1 - C2).flatten() if switch_shape in (5,50,51) \
                        else _depth_parabolic(depth_vec, Hw).flatten()
    cavity_depth_vec = z_cav.copy()  # mid-depths used for cavities

    return horizontal_displacement, transverse_displacement, vertical_displacement, cavity_depth_vec, delta_wall_vector
