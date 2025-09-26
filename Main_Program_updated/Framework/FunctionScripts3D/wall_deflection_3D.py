import numpy as np
from typing import Iterable, Tuple, Union, Optional, Dict

# -----------------------------
#  Kernel: Eq_shaft_3d_AF_v2
# -----------------------------
def eq_shaft_3d(x: np.ndarray,
                y: np.ndarray,
                z: np.ndarray,
                h: float,
                nu: float,
                a: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Port of MATLAB Eq_shaft_3d_AF_v2 (spherical cavity + image terms).

    Parameters
    ----------
    x, y, z : arrays (same shape)
        Field coordinates where displacements are computed.
    h : float
        Depth of the cavity center (z positive downward).
    nu : float
        Poisson's ratio of the soil.
    a : float
        Equivalent spherical cavity radius.

    Returns
    -------
    ux, uy, uz : arrays (same shape as x)
        Displacements, with uz positive downward to match the MATLAB sign convention.
    """
    R1 = np.sqrt(x**2 + y**2 + (z - h)**2)
    R2 = np.sqrt(x**2 + y**2 + (z + h)**2)

    # Avoid division by zero (very close to the cavity center we cap radii)
    eps = np.finfo(float).eps
    R1 = np.maximum(R1, eps)
    R2 = np.maximum(R2, eps)

    common = (1.0 / R1**3) + (3.0 - 4.0*nu) / R2**3 - 6.0*z*(z + h) / R2**5

    ux = - (1.0/3.0) * a**3 * x * common
    uy = - (1.0/3.0) * a**3 * y * common

    fz = ( - (z - h) / R1**3
           - 2.0*z / R2**3
           + (3.0 - 4.0*nu) * (z + h) / R2**3
           + 6.0*z*(z + h)**2 / R2**5 )
    uz = (1.0/3.0) * a**3 * fz
    return ux, uy, uz

# ---------------------------------------------------
#  Vertical wall-deflection shapes (per depth h = z)
# ---------------------------------------------------
def delta_wall_parabola(beta: float, h: np.ndarray, Hw: float) -> np.ndarray:
    """Parabolic shape: -6*(beta/Hw)*h^2 + 6*beta*h (β normalized by Hw)."""
    return -6.0 * beta / Hw * h**2 + 6.0 * beta * h

def delta_wall_cantilever(beta: float, h: np.ndarray, Hw: float) -> np.ndarray:
    """Cantilever shape: 2*beta*Hw - 2*beta*h."""
    return  2.0 * beta * Hw - 2.0 * beta * h

def delta_wall_kickin(beta: float, h: np.ndarray, Hw: float) -> np.ndarray:
    """Kick-in shape (as in your 2D code): -4*(beta/Hw)*h^2 + 4.7*beta*h."""
    return -4.0 * beta / Hw * h**2 + 4.7 * beta * h

def delta_wall_custom_C1C2C3(beta: float,
                             h: np.ndarray,
                             Hw: float,
                             C1: float,
                             C2: float) -> np.ndarray:
    """
    Custom linear combination: C1*cantilever + C2*parabolic + C3*kick-in, with C3 = 1 - C1 - C2.
    (Matches Wall_deflection_v_38.py shape=5 logic.)
    """
    C3 = 1.0 - C1 - C2
    if not np.isclose(C1 + C2 + C3, 1.0, atol=1e-8):
        raise ValueError("(C1 + C2 + C3) must equal 1.0")
    return (C1 * delta_wall_cantilever(beta, h, Hw)
            + C2 * delta_wall_parabola(beta, h, Hw)
            + C3 * delta_wall_kickin(beta, h, Hw))


# ------------------------------------------------------------
#  Longitudinal attenuation along wall (y_midwall coordinate)
# ------------------------------------------------------------
def atten_none(y_midwall: np.ndarray,
               Hw: float,
               He_Hw: float,
               Lwall: float) -> np.ndarray:
    """No attenuation; returns ones with same shape as y_midwall."""
    return np.ones_like(y_midwall, dtype=float)

def atten_mu_huang_2016(y_midwall: np.ndarray,
                        Hw: float,
                        He_Hw: float,
                        Lwall: float) -> np.ndarray:
    """
    Mu & Huang (2016) Gaussian-like attenuation used in your MATLAB (switch_shape==30).
    The MATLAB expression reduces to exp(-pi * ( y / (0.5*Lwall*(0.069*log(He/Lwall)+1.03)) )^2 ),
    where He = He_Hw * Hw.
    """
    He = He_Hw * Hw
    denom = 0.5 * Lwall * (0.069 * np.log(He / Lwall) + 1.03)
    # guard small/negative denom due to logs in edge cases
    denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
    return np.exp(-1.0 * np.pi * (y_midwall / denom)**2)


# ------------------------------------
#  Utility: build rectangle perimeter
# ------------------------------------
def perimeter_sample_points(Lx: float, Ly: float, delta_perimeter: float):
    """
    Walk the rectangle perimeter (bottom -> right -> top -> left) and return:
    - arrays Xc_vec, Yc_vec with equally spaced points;
    - index arrays for each wall;
    - per-wall effective lengths (2*Lx or 2*Ly);
    - straight "y_midwall" coordinate per wall (along-wall coordinate).
    """
    P = 2.0 * (2.0*Lx + 2.0*Ly)
    num_stacks = max(1, int(np.ceil(P / delta_perimeter)))
    s = np.linspace(0.0, P, num_stacks, endpoint=False)

    Xc_vec = np.zeros_like(s)
    Yc_vec = np.zeros_like(s)

    # bottom: y = -Ly, x from -Lx -> +Lx (length = 2Lx)
    idx1 = np.where(s < 2.0*Lx)[0]

    # right: x = +Lx, y from -Ly -> +Ly (length = 2Ly)
    idx2 = np.where((s >= 2.0*Lx) & (s < 2.0*Lx + 2.0*Ly))[0]

    # top: y = +Ly, x from +Lx -> -Lx (length = 2Lx)
    idx3 = np.where((s >= 2.0*Lx + 2.0*Ly) & (s < 4.0*Lx + 2.0*Ly))[0]

    # left: x = -Lx, y from +Ly -> -Ly (length = 2Ly)
    idx4 = np.where(s >= 4.0*Lx + 2.0*Ly)[0]

    # Assign coordinates
    if idx1.size:
        Xc_vec[idx1] = -Lx + s[idx1]
        Yc_vec[idx1] = -Ly
    if idx2.size:
        Xc_vec[idx2] = +Lx
        Yc_vec[idx2] = -Ly + (s[idx2] - 2.0*Lx)
    if idx3.size:
        Xc_vec[idx3] = +Lx - (s[idx3] - (2.0*Lx + 2.0*Ly))
        Yc_vec[idx3] = +Ly
    if idx4.size:
        Xc_vec[idx4] = -Lx
        Yc_vec[idx4] = +Ly - (s[idx4] - (4.0*Lx + 2.0*Ly))

    Lwall = np.array([2.0*Lx, 2.0*Ly, 2.0*Lx, 2.0*Ly], dtype=float)
    return (Xc_vec, Yc_vec,
            idx1, idx2, idx3, idx4,
            Lwall, s)


# ------------------------------------------------------
#  Main driver (XY plane at z=0 by default; surface)
# ------------------------------------------------------
def wall_deflection_3D(
    Hw: float,
    Lx: float,
    Ly: float,
    nu: float = 0.499,
    beta_walls: Union[float, Iterable[float]] = 0.075/100.0,
    vertical_shape: str = "parabola",       # "parabola" | "custom"
    longitudinal: str = "none",             # "none" | "mu_huang_2016"
    He_Hw: float = 0.5,                     # He/Hw ratio (for longitudinal shapes)
    C1: Optional[float] = None,             # for vertical_shape == "custom"
    C2: Optional[float] = None,             # for vertical_shape == "custom"
    delta_z_cavities: Optional[float] = None,
    delta_perimeter: float = 2.5,
    spacing_x: Optional[float] = None,
    spacing_y: Optional[float] = None,
    spacing_z: Optional[float] = None,
    plane: str = "xy",                      # "xy" | "xz" | "xyz" (xyz not yet implemented fully)
    extent_xy: float = 5.0,                 # domain extent multiplier of Hw
    solution_type: int = 2,                 # 1=single wall (m=2), 2=four walls (m=1)
    zero_inside: bool = True
) -> Dict[str, np.ndarray]:
    """
    Compute displacements from a 3D rectangular station by superposition of spherical cavities.

    Parameters match your MATLAB defaults where possible. Start simple and iterate.

    Returns
    -------
    dict with X, Y, Z, ux, uy, uz, and meta.
    """
    # --- meshing defaults
    if delta_z_cavities is None:
        delta_z_cavities = Hw / 10.0
    if spacing_x is None:
        spacing_x = Hw * 0.05
    if spacing_y is None:
        spacing_y = Hw * 0.05
    if spacing_z is None:
        spacing_z = Hw * 0.05

    # --- output grid (start with surface plane z=0 or vertical xz plane)
    if plane == "xy":
        x_vect = np.arange(-extent_xy*Hw, extent_xy*Hw + spacing_x, spacing_x)
        y_vect = np.arange(-extent_xy*Hw, extent_xy*Hw + spacing_y, spacing_y)
        z_vect = np.array([0.0])
        Xsection, Ysection = np.meshgrid(x_vect, y_vect, indexing="xy")
        Zsection = np.zeros_like(Xsection)
    elif plane == "xz":
        x_vect = np.arange(+Lx + spacing_x, extent_xy*Hw + spacing_x, spacing_x)  # as in MATLAB example
        z_vect = np.arange(0.0, 2.0*Hw + spacing_z, spacing_z)
        Xsection, Zsection = np.meshgrid(x_vect, z_vect, indexing="xy")
        Ysection = np.zeros_like(Xsection)
    else:
        raise NotImplementedError("plane='xyz' not implemented yet. Use 'xy' or 'xz'.")

    # --- depth discretization
    z_edges = np.arange(0.0, Hw + delta_z_cavities, delta_z_cavities)
    if z_edges[-1] < Hw - 1e-12:
        z_edges = np.append(z_edges, Hw)
    z_cav_vec = 0.5 * (z_edges[:-1] + z_edges[1:])  # cavity centers along depth

    # --- perimeter sampling
    (Xc_vec, Yc_vec, idx1, idx2, idx3, idx4, Lwall, s) = perimeter_sample_points(Lx, Ly, delta_perimeter)

    # --- per-wall beta (accept scalar or iterable of 4)
    if np.isscalar(beta_walls):
        beta_arr = np.array([beta_walls]*4, dtype=float)
    else:
        arr = np.array(list(beta_walls), dtype=float)
        if arr.size != 4:
            raise ValueError("beta_walls must be a scalar or an iterable of length 4.")
        beta_arr = arr

    # --- vertical shape function
    def vertical_delta(beta, h):
        if vertical_shape == "parabola":
            return delta_wall_parabola(beta, h, Hw)
        elif vertical_shape == "custom":
            if C1 is None or C2 is None:
                raise ValueError("Provide C1 and C2 when vertical_shape='custom'. (C3 = 1 - C1 - C2)")
            return delta_wall_custom_C1C2C3(beta, h, Hw, C1, C2)
        else:
            raise ValueError("vertical_shape must be 'parabola' or 'custom' for now.")

    # --- longitudinal attenuation
    def long_atten(y_midwall, Lwall_eff):
        if longitudinal == "none":
            return atten_none(y_midwall, Hw, He_Hw, Lwall_eff)
        elif longitudinal == "mu_huang_2016":
            return atten_mu_huang_2016(y_midwall, Hw, He_Hw, Lwall_eff)
        else:
            raise ValueError("longitudinal must be 'none' or 'mu_huang_2016' for now.")

    # --- build per-wall y_midwall arrays matching MATLAB grids for DELTA_wall and DELTA_cav
    # For each wall we need a through-depth grid and along-wall coordinates at both z-edges and z-centers
    def wall_y_midwall(idx, along_is_x: bool) -> np.ndarray:
        """Return the along-wall coordinate at the sampling points s[idx]."""
        if idx.size == 0:
            return np.array([])
        if along_is_x:  # bottom/top walls: along x
            return -Lx + (s[idx]) if along_is_x else None  # not used
        else:          # right/left walls: along y
            return -Ly + (s[idx] - 2.0*Lx)

    # 2D grids for wall z-edges and cavity z-centers
    # Bottom (1) and Top (3): y fixed, along-wall coordinate is x
    # Right (2) and Left (4): x fixed, along-wall coordinate is y
    Ymid_1_edges = -Lx + s[idx1] if idx1.size else np.array([])
    Ymid_3_edges = +Lx - (s[idx3] - (2.0*Lx + 2.0*Ly)) if idx3.size else np.array([])
    Ymid_2_edges = -Ly + (s[idx2] - 2.0*Lx) if idx2.size else np.array([])
    Ymid_4_edges = +Ly - (s[idx4] - (4.0*Lx + 2.0*Ly)) if idx4.size else np.array([])

    # Create 2D grids (depth x along-wall)
    def mesh_if_any(y_edges):
        if y_edges.size == 0:
            return np.empty((0,0)), np.empty((0,0))
        return np.meshgrid(y_edges, z_edges, indexing="xy")

    def mesh_if_any_centers(y_edges):
        if y_edges.size == 0:
            return np.empty((0,0)), np.empty((0,0))
        return np.meshgrid(y_edges, z_cav_vec, indexing="xy")

    Ymid_wall_1, Z_wall_edges_1 = mesh_if_any(Ymid_1_edges)
    Ymid_wall_2, Z_wall_edges_2 = mesh_if_any(Ymid_2_edges)
    Ymid_wall_3, Z_wall_edges_3 = mesh_if_any(Ymid_3_edges)
    Ymid_wall_4, Z_wall_edges_4 = mesh_if_any(Ymid_4_edges)

    Ymid_cav_1, Z_cav_1 = mesh_if_any_centers(Ymid_1_edges)
    Ymid_cav_2, Z_cav_2 = mesh_if_any_centers(Ymid_2_edges)
    Ymid_cav_3, Z_cav_3 = mesh_if_any_centers(Ymid_3_edges)
    Ymid_cav_4, Z_cav_4 = mesh_if_any_centers(Ymid_4_edges)

    # --- compute DELTA (wall deflection) at z-edges (for volume per run-meter) and at z-centers (for cavities)
    def wall_delta_grid(beta_val, Ymid_edges, Z_edges, Ymid_centers, Z_centers, Lwall_eff):
        if Ymid_edges.size == 0:
            return (np.empty((0,0)), np.empty((0,0)))
        # vertical component
        delta_edges = vertical_delta(beta_val, Z_edges, )
        delta_centers = vertical_delta(beta_val, Z_centers, )
        # longitudinal attenuation
        atten_edges = long_atten(Ymid_edges, Lwall_eff)
        atten_centers = long_atten(Ymid_centers, Lwall_eff)
        return (delta_edges * atten_edges, delta_centers * atten_centers)

    # Fix vertical_delta signature: expects (beta, h)
    def vertical_delta(beta_val, h_vals):
        return (delta_wall_parabola(beta_val, h_vals, Hw) if vertical_shape == "parabola"
                else delta_wall_custom_C1C2C3(beta_val, h_vals, Hw, C1, C2))

    DELTA_wall_1, DELTA_cav_1 = wall_delta_grid(beta_arr[0], Ymid_wall_1, Z_wall_edges_1, Ymid_cav_1, Z_cav_1, Lwall[0])
    DELTA_wall_2, DELTA_cav_2 = wall_delta_grid(beta_arr[1], Ymid_wall_2, Z_wall_edges_2, Ymid_cav_2, Z_cav_2, Lwall[1])
    DELTA_wall_3, DELTA_cav_3 = wall_delta_grid(beta_arr[2], Ymid_wall_3, Z_wall_edges_3, Ymid_cav_3, Z_cav_3, Lwall[2])
    DELTA_wall_4, DELTA_cav_4 = wall_delta_grid(beta_arr[3], Ymid_wall_4, Z_wall_edges_4, Ymid_cav_4, Z_cav_4, Lwall[3])

    # --- per-run-meter volume loss (for debugging/consistency)
    def trapz_depth(delta_edges):
        """
        Integrate along the depth dimension using z_edges as the x grid.
        We auto-detect which axis matches len(z_edges) to avoid shape mismatches.
        """
        if delta_edges.size == 0:
            return np.array([])
        return np.array([])
        # pick the axis whose length matches z_edges
        if delta_edges.shape[-1] == z_edges.size:
            axis_depth = -1
        elif delta_edges.shape[0] == z_edges.size:
            axis_depth = 0
        else:
            # fall back: assume depth is the second axis for (along, depth) grids
            axis_depth = 1 if delta_edges.ndim > 1 else 0
        return np.trapz(delta_edges, x=z_edges, axis=axis_depth)



    VLW_1 = trapz_depth(DELTA_wall_1)
    VLW_2 = trapz_depth(DELTA_wall_2)
    VLW_3 = trapz_depth(DELTA_wall_3)
    VLW_4 = trapz_depth(DELTA_wall_4)

    beta_runmeter_1 = VLW_1 / Hw**2 if VLW_1.size else np.array([])
    beta_runmeter_2 = VLW_2 / Hw**2 if VLW_2.size else np.array([])
    beta_runmeter_3 = VLW_3 / Hw**2 if VLW_3.size else np.array([])
    beta_runmeter_4 = VLW_4 / Hw**2 if VLW_4.size else np.array([])

    # --- superposition loop
    ux = np.zeros_like(Xsection, dtype=float)
    uy = np.zeros_like(Ysection, dtype=float)
    uz = np.zeros_like(Zsection, dtype=float)

    # symmetry multiplier
    if solution_type == 1:
        m_sym = 2.0  # single wall
    elif solution_type == 2:
        m_sym = 1.0  # four walls analytical
    else:
        raise NotImplementedError("solution_type=3 (semi-analytical taper) to be added next.")

    def ensure_along_depth(arr: np.ndarray, depth_len: int) -> np.ndarray:
        """
        Ensure array is shaped (along, depth). If it's (depth, along),
        transpose it. If it's empty, return as-is.
        """
        if arr.size == 0 or arr.ndim != 2:
            return arr
        a0, a1 = arr.shape
        if a1 == depth_len:
            return arr                   # already (along, depth)
        if a0 == depth_len:
            return arr.T                 # was (depth, along) -> transpose
        # Fallback: if neither axis matches, keep as-is (we also guard in the loop)
        return arr





    # helper to iterate over all stacks in a wall
    def add_wall_contrib(idx, Xc, Yc, DELTA_cav_wall, which_wall, Lwall_eff):
        nonlocal ux, uy, uz
        if idx.size == 0 or DELTA_cav_wall.size == 0:
            return

        for j_local, j in enumerate(idx):
            deltaXl = Xc[j]
            deltaYl = Yc[j]

            # Depth profile for this along-wall stack
            if DELTA_cav_wall.ndim == 2 and j_local < DELTA_cav_wall.shape[0]:
                delta_col = DELTA_cav_wall[j_local, :]  # shape (depth,)
            else:
                continue  # nothing to add for this stack

            # radius per depth slice
            a_vec = (m_sym * 3.0/(4.0*np.pi) * delta_col * delta_z_cavities * delta_perimeter) ** (1.0/3.0)

            # Loop over the shorter of (len(z_cav_vec), len(a_vec))
            depth_len = min(z_cav_vec.size, a_vec.size)
            if depth_len == 0:
                continue

            Xloc = Xsection - deltaXl
            Yloc = Ysection - deltaYl
            Zloc = Zsection

            for k in range(depth_len):
                a_k = float(a_vec[k])
                if not np.isfinite(a_k) or a_k <= 0.0:
                    continue
                h = float(z_cav_vec[k])
                uxi, uyi, uzi = eq_shaft_3d(Xloc, Yloc, Zloc, h, nu, a_k)
                ux += uxi
                uy += uyi
                uz += uzi

    add_wall_contrib(idx1, Xc_vec, Yc_vec, DELTA_cav_1 if DELTA_cav_1.size else DELTA_cav_1, which_wall=1, Lwall_eff=Lwall[0])
    add_wall_contrib(idx2, Xc_vec, Yc_vec, DELTA_cav_2 if DELTA_cav_2.size else DELTA_cav_2, which_wall=2, Lwall_eff=Lwall[1])
    add_wall_contrib(idx3, Xc_vec, Yc_vec, DELTA_cav_3 if DELTA_cav_3.size else DELTA_cav_3, which_wall=3, Lwall_eff=Lwall[2])
    add_wall_contrib(idx4, Xc_vec, Yc_vec, DELTA_cav_4 if DELTA_cav_4.size else DELTA_cav_4, which_wall=4, Lwall_eff=Lwall[3])

    # --- zero out inside the station footprint (optional)
    if zero_inside:
        inside = (Xsection >= -Lx) & (Xsection <= Lx) & (Ysection >= -Ly) & (Ysection <= Ly)
        ux = np.where(inside, 0.0, ux)
        uy = np.where(inside, 0.0, uy)
        uz = np.where(inside, 0.0, uz)

    meta = dict(
        z_edges=z_edges, z_cav_vec=z_cav_vec, perimeter_samples=len(s),
        VLW_per_run_meter=[beta_runmeter_1, beta_runmeter_2, beta_runmeter_3, beta_runmeter_4],
        Lwall=Lwall, Xc_vec=Xc_vec, Yc_vec=Yc_vec,
        beta_walls=beta_arr, vertical_shape=vertical_shape, longitudinal=longitudinal,
        solution_type=solution_type, He_Hw=He_Hw,
        delta_z_cavities=delta_z_cavities, delta_perimeter=delta_perimeter
    )

    return dict(X=Xsection, Y=Ysection, Z=Zsection, ux=ux, uy=uy, uz=uz, meta=meta)


# --------------------------
#  Tiny self-test (optional)
# --------------------------
if __name__ == "__main__":
    # A very small smoke test to ensure the function runs without error.
    out = wall_deflection_3D(
        Hw=19.0,
        Lx=9.5,
        Ly=32.0,
        nu=0.499,
        beta_walls=[0.075/100]*4,
        vertical_shape="parabola",
        longitudinal="none",
        He_Hw=0.5,
        delta_z_cavities=19.0/10.0,
        delta_perimeter=2.5,
        spacing_x=19.0*0.2,
        spacing_y=19.0*0.2,
        plane="xy",
        extent_xy=1.5,
        solution_type=2,
        zero_inside=True
    )
    # Print a quick summary
    ux, uy, uz = out["ux"], out["uy"], out["uz"]
    print("Grid shape:", ux.shape)
    print("Max |ux|, |uy|, |uz| [m]:", float(np.nanmax(np.abs(ux))), float(np.nanmax(np.abs(uy))), float(np.nanmax(np.abs(uz))))
    print("Done.")
