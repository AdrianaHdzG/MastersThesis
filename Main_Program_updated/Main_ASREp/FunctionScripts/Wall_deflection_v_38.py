import numpy as np


# def wall_deflection(retaining_wall_depth, avg_wall_displacement, wall_deflection_shape, soil_ovalization, volumetric, x_coords=None):
def wall_deflection(retaining_wall_depth, avg_wall_displacement, wall_deflection_shape, soil_ovalization, volumetric,
                    x_coords=None, C1_val=None, C2_val=None, foundation_depth=None):
    """
    Calculates the wall deflection based on the given parameters.

    Parameters:
    retaining_wall_depth : float       - Depth of the retaining wall [m]
    avg_wall_displacement : float      - Average wall displacement normalized by the wall height [-]
    wall_deflection_shape : int        - Shape of wall deflection M0-M4 [-]
    soil_ovalization : float           - Soil ovalization [-]
    volumetric : float                 - Volumetric (1 = no contraction and expansion) [-]
    x_coords : array-like, optional    - Custom x-coordinates for deflection calculation.
                                         If provided, overrides the default x_coords.
    C1_val
    C2_val

    Returns:
    horizontal_displacement, vertical_displacement, wall_deflection_all, volume_loss_wall,
    max_wall_deflection, check_ratio, volume_loss_surface, ratio_volume
    """
    # print(retaining_wall_depth, avg_wall_displacement, wall_deflection_shape, soil_ovalization, volumetric) DEBUGGING
    deltaH = retaining_wall_depth / 100  # Integration step height

    # Define x_coords if not provided
    if x_coords is None:
        grid_perH = 0.01
        x_coords = np.arange(retaining_wall_depth * 0.01, 2.0 * retaining_wall_depth, retaining_wall_depth * grid_perH)

    # Define the wall deflection functions based on the shape
    def wall_uniform_deflection(depth):
        return 1 * avg_wall_displacement * retaining_wall_depth + 0 * depth

    def wall_cantilever_deflection(depth):
        return 2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth

    def wall_parabola_deflection(depth):
        return -6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth

    def wall_composite_deflection(depth):
        return 0.35 * (2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth) + \
            0.65 * (-6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth)

    def wall_kickin_deflection(depth):
        return -4 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 4.7 * avg_wall_displacement * depth

    def wall_custom_deflection(depth):  # C1, C2, C5

        C3_val = 1 - C1_val - C2_val
        if (C1_val + C2_val + C3_val) < 0.99999 or (C1_val + C2_val + C3_val) > 1.00001:
            print('(C1 + C2 + C3) != 1, this is not consistent with assumptions')
            raise ValueError('(C1 + C2 + C3) != 1, this is not consistent with assumptions', C1_val, C2_val, C3_val)
        return C1_val * (2 * avg_wall_displacement * retaining_wall_depth - 2 * avg_wall_displacement * depth) + \
            C2_val * (
                        -6 * avg_wall_displacement / retaining_wall_depth * depth ** 2 + 6 * avg_wall_displacement * depth) + \
            C3_val * (2 * avg_wall_displacement * depth)

    # Select the deflection function based on wall_deflection_shape
    if wall_deflection_shape == 0:
        wall_deflection_func = wall_uniform_deflection
    elif wall_deflection_shape == 1:
        wall_deflection_func = wall_cantilever_deflection
    elif wall_deflection_shape == 2:
        wall_deflection_func = wall_parabola_deflection
    elif wall_deflection_shape == 3:
        wall_deflection_func = wall_composite_deflection
    elif wall_deflection_shape == 4:
        wall_deflection_func = wall_kickin_deflection
    elif wall_deflection_shape == 5:
        wall_deflection_func = wall_custom_deflection
    else:
        if wall_deflection_shape is not None and not isinstance(wall_deflection_shape, (int, np.ndarray)):
            raise ValueError("Error: wall_deflection_shape must be an integer or numpy array")
        else:
            raise ValueError("Error: wall_deflection_shape must be either 0, 1, 2, 3, 4, 5")

    if foundation_depth is None:
        z_coords = np.array([0])  # Only need to calculate surface, so this can be set to 0
    else:
        z_coords = np.array([foundation_depth])  # Only need to calculate a given depth

    depth_vec = np.arange(0, retaining_wall_depth + deltaH, deltaH)  # Integrates over wall depth, could go from found_d
    cavity_depth_vec = np.arange(deltaH / 2, retaining_wall_depth - deltaH / 2 + deltaH, deltaH)
    num_cavities = len(cavity_depth_vec)

    # Create the mesh grid
    x_section, z_section = np.meshgrid(x_coords, z_coords)

    # Superposition method (we don't need this for loop)
    horizontal_displacement = np.zeros_like(x_section)
    vertical_displacement = np.zeros_like(x_section)
    deflection_vec = wall_deflection_func(cavity_depth_vec)
    sign_vector = np.sign(deflection_vec)
    radius_vector = np.sqrt(np.abs(deflection_vec * deltaH / np.pi))  # Refer to GF paper or thesis

    # Calculate displacement increments and superposition
    for i in range(num_cavities):
        epsilon = sign_vector[i]
        radius = radius_vector[i]
        cavity_depth = cavity_depth_vec[i]
        u_xinc, u_zinc = u_GON(z_section, x_section, cavity_depth, radius, epsilon, soil_ovalization * epsilon,
                               volumetric)  # u_GON(z, x, depth, radius, epsilon, soil_ovalization, volumetric)
        horizontal_displacement += u_xinc
        vertical_displacement += u_zinc

    wall_deflection_all = wall_deflection_func(depth_vec)
    volume_loss_wall = np.trapz(wall_deflection_all, depth_vec)  # Volume loss of wall deflection
    max_wall_deflection = np.max(wall_deflection_all)  # Maximum wall deflection
    volume_loss_wall_v2 = avg_wall_displacement * retaining_wall_depth ** 2
    check_ratio = volume_loss_wall / volume_loss_wall_v2

    volume_loss_surface = np.trapz(vertical_displacement[0, :], x_section[0, :])  # Volume loss of wall deflection
    ratio_volume = volume_loss_surface / volume_loss_wall

    return horizontal_displacement, vertical_displacement, wall_deflection_all, volume_loss_wall, max_wall_deflection, check_ratio, volume_loss_surface, ratio_volume


def u_GON(z, x, depth, radius, epsilon, soil_ovalization, volumetric):
    """
    Calculates the displacement increments.

    Parameters:
    z : array_like        - Spatial coordinate z
    x : array_like        - Spatial coordinate x
    depth : float         - Tunnel depth
    radius : float        - Tunnel radius
    epsilon : float       - Equal to Vlt./2
    soil_ovalization : float   - Parameter for ovalization
    volumetric : float    - Equal to 1 for elastic behaviour

    Returns:
    ux : array_like       - Displacement in x direction
    uz : array_like       - Displacement in z direction
    """

    ovalization_ratio = soil_ovalization / epsilon
    normalized_x = x / depth
    depth_below = z - depth
    depth_above = z + depth
    normalized_depth_below = depth_below / depth
    normalized_depth_above = depth_above / depth
    normalized_depth = z / depth
    radius_below = np.sqrt(x ** 2 + depth_below ** 2) / depth
    radius_above = np.sqrt(x ** 2 + depth_above ** 2) / depth

    # General expression
    ux = 2 * epsilon * radius * (radius / depth) ** (2 * volumetric - 1) * (
            -((normalized_x * (1 - (
                    ovalization_ratio * (normalized_x ** 2 - normalized_depth_below ** 2)) / radius_below ** 2)) / (
                      2 * radius_below ** (2 * volumetric))) -
            (normalized_x * (1 - (
                    ovalization_ratio * (normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2)) / (
                    2 * radius_above ** (2 * volumetric)) +
            (4 * normalized_x * normalized_depth * (normalized_depth_above / radius_above ** 2 - (ovalization_ratio * (
                    normalized_x ** 2 - 3 * normalized_depth_above ** 2)) / radius_above ** 4)) / (
                    2 * radius_above ** (2 * volumetric)))

    uz = 2 * epsilon * radius * (radius / depth) ** (2 * volumetric - 1) * (
            -((normalized_depth_below * (1 - (
                    ovalization_ratio * (normalized_x ** 2 - normalized_depth_below ** 2)) / radius_below ** 2)) / (
                      2 * radius_below ** (2 * volumetric))) +
            (normalized_depth_above * (1 + (
                    ovalization_ratio * (normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2)) / (
                    2 * radius_above ** (2 * volumetric)) -
            ((2 * (normalized_depth + ovalization_ratio) * (
                    normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 2 +
             (4 * ovalization_ratio * normalized_depth * normalized_depth_above * (
                     3 * normalized_x ** 2 - normalized_depth_above ** 2)) / radius_above ** 4) / (
                    2 * radius_above ** (2 * volumetric)))

    return ux, uz


# Example call to the function with the initial parameters
# horizontal_displacement, vertical_displacement, wall_deflection_all, volume_loss_wall, max_wall_deflection, check_ratio, volume_loss_surface, ratio_volume = wall_deflection(6, 0.01, 2, 0, 1)


def wall_deflection_direct(retaining_wall_depth, soil_ovalization, volumetric, delta_wall, depth_w,
                           x_coords=None, foundation_depth=None):
    """
    Calculates the wall deflection based on the given parameters.

    Parameters:
    retaining_wall_depth : float       - Depth of the retaining wall [m]
    avg_wall_displacement : float      - Average wall displacement normalized by the wall height [-]
    wall_deflection_shape : int        - Shape of wall deflection M0-M4 [-]
    soil_ovalization : float           - Soil ovalization [-]
    volumetric : float                 - Volumetric (1 = no contraction and expansion) [-]
    x_coords : array-like, optional    - Custom x-coordinates for deflection calculation.
                                         If provided, overrides the default x_coords.
    C1_val
    C2_val

    Returns:
    horizontal_displacement, vertical_displacement, wall_deflection_all, volume_loss_wall,
    max_wall_deflection, check_ratio, volume_loss_surface, ratio_volume
    """
    # print(retaining_wall_depth, avg_wall_displacement, wall_deflection_shape, soil_ovalization, volumetric) DEBUGGING
    deltaH = retaining_wall_depth / 100  # Integration step height

    # Define x_coords if not provided
    if x_coords is None:
        grid_perH = 0.01
        x_coords = np.arange(retaining_wall_depth * 0.01, 2.0 * retaining_wall_depth, retaining_wall_depth * grid_perH)

    if foundation_depth is None:
        z_coords = np.array([0])  # Only need to calculate surface, so this can be set to 0
    else:
        z_coords = np.array([foundation_depth])  # Only need to calculate a given depth

    depth_vec = np.arange(0, retaining_wall_depth + deltaH, deltaH)  # Integrates over wall depth, could go from found_d
    cavity_depth_vec = np.arange(deltaH / 2, retaining_wall_depth - deltaH / 2 + deltaH, deltaH)
    num_cavities = len(cavity_depth_vec)

    # Interpolate the data
    from scipy.interpolate import PchipInterpolator
    # Create the interpolator
    pchip_interp = PchipInterpolator(depth_w, delta_wall)
    # Perform interpolation
    delta_wall_vector = pchip_interp(cavity_depth_vec)

    # Create the mesh grid
    x_section, z_section = np.meshgrid(x_coords, z_coords)

    # Superposition method (we don't need this for loop)
    horizontal_displacement = np.zeros_like(x_section)
    vertical_displacement = np.zeros_like(x_section)
    # deflection_vec = wall_deflection_func(cavity_depth_vec)
    sign_vector = np.sign(delta_wall_vector)
    radius_vector = np.sqrt(np.abs(delta_wall_vector * deltaH / np.pi))  # Refer to GF paper or thesis

    # Calculate displacement increments and superposition
    for i in range(num_cavities):
        epsilon = sign_vector[i]
        radius = radius_vector[i]
        cavity_depth = cavity_depth_vec[i]
        u_xinc, u_zinc = u_GON(z_section, x_section, cavity_depth, radius, epsilon, soil_ovalization * epsilon,
                               volumetric)  # u_GON(z, x, depth, radius, epsilon, soil_ovalization, volumetric)
        horizontal_displacement += u_xinc
        vertical_displacement += u_zinc

    # wall_deflection_all = wall_deflection_func(depth_vec)
    # volume_loss_wall = np.trapz(wall_deflection_all, depth_vec)  # Volume loss of wall deflection
    # max_wall_deflection = np.max(wall_deflection_all)  # Maximum wall deflection
    # volume_loss_wall_v2 = avg_wall_displacement * retaining_wall_depth ** 2
    # check_ratio = volume_loss_wall / volume_loss_wall_v2

    # volume_loss_surface = np.trapz(vertical_displacement[0, :], x_section[0, :])  # Volume loss of wall deflection
    # ratio_volume = volume_loss_surface / volume_loss_wall

    return horizontal_displacement, vertical_displacement, cavity_depth_vec, delta_wall_vector
