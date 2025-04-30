import numpy as np


def categorize_damage(element_strain):
    """
    Categorize the damage based on the computed tensile strain.

    Parameters:
    tensile_strain : float - Computed tensile strain

    Returns:
    damage_category : str - The damage category based on the tensile strain
    """

    tensile_strain = element_strain * 100  # Convert to percent
    print("tensile strain is %", tensile_strain)
    if tensile_strain < 0.05:
        return 0, tensile_strain
    elif tensile_strain < 0.075:
        return 1, tensile_strain
    elif tensile_strain < 0.15:
        return 2, tensile_strain
    elif tensile_strain < 0.3:
        return 3, tensile_strain
    elif tensile_strain < 0.6:
        return 4, tensile_strain
    else:
        return 5, tensile_strain


def compute_tensile_strain_Jinyan(model, model_properties):
    """
    Compute the tensile strain for a beam element.

    Parameters:
    axial_elongation : float - Axial elongation of the beam element
    bending_curvature : float - Bending curvature of the beam element
    element_length : float - Length of the beam element

    BASED ON AGE_lecture_3_4 from Andrea Franza

    Returns:
    tensile_strain : float - Computed tensile strain
    """

    EA = model_properties[0]
    EI = model_properties[1]
    Gb = model_properties[2]
    Ab = model_properties[3]
    poissons_ratio = model_properties[4]
    shear_factor_midpoint = model_properties[5]
    building_height = model_properties[6]
    # dist_NA = model_properties[7]
    neutral_line = model_properties[8]

    As = Ab * (10 + (poissons_ratio * 10)) / (12 + (11 * poissons_ratio))  # According to Wikipedia
    GAs = Gb * As

    # N / EA
    strain_axial_normal = model.axialForce / EA

    # Current validation step

    # ( M / EI ) * d
    strain_axial_bending_top = -(model.moment / EI) * (building_height - neutral_line)  # Compressive stains in the top
    strain_axial_bending_bottom = (model.moment / EI) * neutral_line

    # Navier's formula
    tensile_strain_top = strain_axial_normal + strain_axial_bending_top
    tensile_strain_bottom = strain_axial_normal + strain_axial_bending_bottom

    # gamma = V/GAs & e_xy = gamma / 2
    true_shear_strain = (model.shearForce * shear_factor_midpoint / GAs) / 2

    # Bending's contribution e_xx * (1 - v) / 2
    strain_exx_contribution = ((strain_axial_normal * (1 - poissons_ratio)) / 2)
    # Formula 7 from slides
    tensile_strain_midpoint = strain_exx_contribution + np.sqrt(
        (((strain_axial_normal * (1 + poissons_ratio)) / 2) ** 2)
        + (true_shear_strain ** 2))

    dataReturn = {  # Make a DATA struct for all models run
        'axial_strain_exx': strain_axial_normal,
        'max_exx': max(strain_axial_normal),
        'strain_axial_bending_top_exx,b': strain_axial_bending_top,
        'strain_axial_bending_bottom_exx,b': strain_axial_bending_bottom,
        'exx,t,b': tensile_strain_top,
        'exx,b,b': tensile_strain_bottom,
        'e_t_b': max(max(tensile_strain_top), max(tensile_strain_bottom)),  # Only max positive strains
        'true_shear_strain': true_shear_strain,
        'tensile_strain_midpoint': tensile_strain_midpoint,
        'max_e_xy': max(true_shear_strain),
        'max_e_t_mid': max(tensile_strain_midpoint),
        'max_e_t': max(max(max(tensile_strain_top), max(tensile_strain_bottom)), max(tensile_strain_midpoint))
    }

    return dataReturn


# %% Function to generate output - should generally not be touched!
def generate_output(dataReturn, max_eps_t_ssi, output_fields, error_flag=0):
    '''
    This function generates a string based on specified output fields and writes it to a file named `results.out`.
    It can handle an error flag to determine whether to return real values or zeros with the last value set to 1.

    Inputs:
    - `dataReturn`: A dictionary containing various data fields and their values.
    - `max_eps_t_ssi`: A list containing the maximum epsilon value.
    - `output_fields`: A list of strings specifying the fields to include in the output.
    - `error_flag` (optional): An integer (default is 0). If set to 0, the function returns real values. If set to 1, the function returns zeros except for the last value, which is set to 1.

    Outputs:
    - `results.out`: A file containing the generated output string based on the specified fields and error flag.
    '''
    # Generate the f.write string based on the specified output fields
    if error_flag == 0:
        output_values = [
            dataReturn[field] if field in dataReturn else max_eps_t_ssi[0] if field == 'max_eps_t_ssi' else 0 for field
            in output_fields]
    else:
        output_values = [0 for _ in output_fields]
        output_values[-1] = 1

    output_string = " ".join(map(str, output_values))

    # Write the output to the results.out file
    with open("results.out", "w") as f:
        f.write(output_string)

