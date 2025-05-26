import numpy as np


# import matplotlib.pyplot as plt
# import os

def strain_analysis_greenfield(disp_vertical, disp_horizontal, length_beam_element, length_beam, building_height,
                               neutral_line, nu, flag, mode):
    """
    Analyze the strain and deformation parameters for greenfield displacements.

    Parameters:
    disp_vertical : array_like - Vertical displacements of the ground [m]
    disp_horizontal : array_like - Horizontal displacements of the ground [m]
    length_beam_element : float - Length of beam element [m]
    length_beam : float - Total length of the beam [m]
    building_height : float - Height of building [m]
    nu : float - Poisson's ratio [-]

    Returns:
    dataReturn
    """
    # save_input_variables(disp_vertical, disp_horizontal, length_beam_element, length_beam, nu)

    if flag == 'WALL':
        # Calculate horizontal strain
        horizontal_strain = np.gradient(disp_horizontal[0], length_beam_element)  # [-] Horizontal strain

        # Calculate tilt
        vert_for_tilt = disp_vertical[0]  # Unpacking, sometimes needed
        omega_greenfield = (vert_for_tilt[-1] - vert_for_tilt[0]) / length_beam  # [rad] tilt of structure

        # Calculate slope
        slope = np.gradient(disp_vertical[0], length_beam_element)
    else:
        # Calculate horizontal strain
        horizontal_strain = np.gradient(disp_horizontal, length_beam_element)  # [-] Horizontal strain

        # Calculate tilt
        omega_greenfield = (disp_vertical[-1] - disp_vertical[0]) / length_beam  # [rad] tilt of structure

        # Calculate slope
        slope = np.gradient(disp_vertical, length_beam_element)

    # %% Calculate angular distortion and diagonal strains
    angular_distortion = abs(slope - omega_greenfield)  # [rad] Angular distortion

    eps_dt = ((horizontal_strain * (1 - nu)) / 2) + np.sqrt((((horizontal_strain * (1 + nu)) / 2) ** 2) +
                                                            ((angular_distortion / 2) ** 2))

    eps_dt_no_tilt = ((horizontal_strain * (1 - nu)) / 2) + np.sqrt((((horizontal_strain * (1 + nu)) / 2) ** 2) +
                                                                    ((abs(slope) / 2) ** 2))

    # %% New section
    if mode == 'Bending':
        # Calculate bending strains based on curvature
        curvature = np.gradient(slope, length_beam_element)
        eps_bending = np.zeros_like(curvature)

        for i in range(len(curvature)):
            if curvature[i] <= 0:  # If curvature is negative (Hogging)

                # Assume that the neutral axis is at building extreme fibre
                eps_bending[i] = -curvature[i] * building_height + horizontal_strain[i]
            else:  # Building is in sagging mode:

                # Neutral line is the distance from beam base to center of bending
                eps_bending[i] = curvature[i] * neutral_line + horizontal_strain[i]

        eps_tensile = np.zeros_like(curvature)
        for i in range(len(eps_tensile)):
            eps_tensile[i] = max(eps_bending[i], eps_dt[i])

        dataReturn = {  # Make a DATA struct for all models run
            'eps_dt': eps_dt,  # shearing strain along the building, length (num_elem) unit [-]
            'eps_dt_max': max(eps_dt),  # max shearing strain along the building, length (1) unit [-]

             'eps_dt_no_tilt': eps_dt_no_tilt,  # Check influence of tilt

            'eps_bt': eps_bending,  # bending strain along the building, length (num_elem) unit [-]
            'eps_bt_max': max(eps_bending),  # bending Tensile strain along the building, length (1) unit [-]

            'eps_t': eps_tensile,  # Tensile strain along the building, length (num_elem) unit [-]
            'eps_t_max': max(eps_tensile),  # max Tensile strain along the building, length (1) unit [-]

            'beta_d': angular_distortion,  # Angular distorsion, length (num_elem) unit [-]
            'beta_d_max': max(angular_distortion),  # max Angular distorsion, length (1) unit [-]

            'eps_h': horizontal_strain,  # Horizontal strain, length (num_elem) unit [-]
            'eps_h_max': max(horizontal_strain),  # max Horizontal strain, length (1) unit [-]

            'S': slope,  # Greenfield slope (uz), length (num_elem) unit [-]
            'kappa': curvature,  # Greenfield curvature unit [1/m]
            'omega_gf': omega_greenfield,  # Tilt of building, length (1) unit [-]
        }
    else:
        dataReturn = {  # Make a DATA struct for all models run
            'eps_dt': eps_dt,  # shearing strain along the building, length (num_elem) unit [-]
            'eps_dt_max': max(eps_dt),  # max shearing strain along the building, length (1) unit [-]

            'eps_dt_no_tilt': eps_dt_no_tilt,  # Check influence of tilt

            'eps_t': eps_dt,  # Tensile strain along the building, length (num_elem) unit [-]
            'eps_t_max': max(eps_dt),  # max Tensile strain along the building, length (1) unit [-]

            'beta_d': angular_distortion,  # Angular distorsion, length (num_elem) unit [-]
            'beta_d_max': max(angular_distortion),  # max Angular distorsion, length (1) unit [-]

            'eps_h': horizontal_strain,  # Horizontal strain, length (num_elem) unit [-]
            'eps_h_max': max(horizontal_strain),  # max Horizontal strain, length (1) unit [-]

            'S': slope,  # Greenfield slope (uz), length (num_elem) unit [-]
            'kappa': curvature,  # Greenfield curvature unit [1/m]
            'omega_gf': omega_greenfield,  # Tilt of building, length (1) unit [-]
        }

    '''
    x = np.linspace(0, length_beam, int(length_beam / length_beam_element + 1))
    plt.plot(x, eps_tensile * 100, 'k-', label='ε_t', linewidth=3)  # Black solid line
    plt.plot(x, horizontal_strain * 100, 'k:', label='ε_h', linewidth=3)  # Black dotted line
    plt.plot(x, (angular_distortion / 2) * 100, 'k--', label='ε_xy', linewidth=3)  # Black dashed line
    plt.plot([min(x), max(x)], [0.0, 0.0], 'k--', label='Cat 0')
    plt.plot([min(x), max(x)], [0.05, 0.05], 'k--', label='Cat 1 Very slight')
    plt.plot([min(x), max(x)], [0.075, 0.075], 'k:', label='Cat 2 Slight')
    plt.plot([min(x), max(x)], [0.15, 0.15], 'k-', label='Cat 3 Moderate')
    plt.ylabel('Deformation (%)')
    plt.xlim([0, length_beam])

    # Move the legend a bit to the right and make the text smaller
    plt.legend(loc='upper left', bbox_to_anchor=(0.64, 1), borderaxespad=0., prop={'size': 10})

    plt.title('Tensile strain along the building due to Greenfield movements')
    plt.xlabel('Offset (m)')

    # Adjust figure size to accommodate the legend within the plot
    plt.gcf().set_size_inches(15 / 2.54, 8 / 2.54)  # Convert cm to inches for figure size
    plt.gcf().set_facecolor('w')
    plt.gca().tick_params(length=0.0200 * 10, width=0.050 * 10)  # Adjust tick length and width

    # Create subfolder with the name in the string inside 'toMATLAB' directory
    save_name = 'Greenfield_standard'
    save_dir = os.path.join('toMATLAB', save_name)
    os.makedirs(save_dir, exist_ok=True)

    # Save the figure as PDF in the subfolder
    plt.savefig(os.path.join(save_dir, save_name + '.pdf'), dpi=300)
    plt.close()
    # Save the data to a txt file in the subfolder
    data_file_path = os.path.join(save_dir, save_name + '_data.txt')
    with open(data_file_path, 'w') as f:
        f.write('x: {}\n'.format(x.tolist()))
        f.write('eps_tensile: {}\n'.format((eps_tensile * 100).tolist()))
        f.write('horizontal_strain: {}\n'.format((horizontal_strain * 100).tolist()))
        f.write('angular_distortion: {}\n'.format((angular_distortion / 2 * 100).tolist()))
    '''

    return dataReturn


def save_input_variables(disp_vertical, disp_horizontal, length_beam_element, length_beam, nu,
                         filename='input_variables.txt'):
    """
    Save the input variables to a file.

    Parameters:
    disp_vertical : array_like - Vertical displacements of the ground [m]
    disp_horizontal : array_like - Horizontal displacements of the ground [m]
    length_beam_element : float - Length of beam element [m]
    length_beam : float - Total length of the beam [m]
    nu : float - Poisson's ratio [-]
    filename : str - Name of the file to save the variables to
    """
    with open(filename, 'w') as f:
        f.write(f"disp_vertical: {disp_vertical.tolist()}\n")
        f.write(f"disp_horizontal: {disp_horizontal.tolist()}\n")
        f.write(f"length_beam_element: {length_beam_element}\n")
        f.write(f"length_beam: {length_beam}\n")
        f.write(f"nu: {nu}\n")

    print(f"Input variables have been saved to {filename}.")


def print_min_max(array):
    print('minimum value')
    print(min(array))
    print('maximum value')
    print(max(array))
