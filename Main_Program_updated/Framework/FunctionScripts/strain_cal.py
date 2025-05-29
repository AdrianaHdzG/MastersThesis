import numpy as np

def compute_tensile_strain(model, model_properties):


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
    dist_NA = model_properties[7]

    As = Ab * (10 + (poissons_ratio * 10)) / (12 + (11 * poissons_ratio))  # According to Wikipedia
    GAs = Gb * As

    # N / EA
    strain_axial_normal = model.axialForce / EA

    # Current validation step

    # ( M / EI ) * d
    strain_axial_bending_top = -(model.moment / EI) * (building_height / 2)  # Compressive stains in the top
    strain_axial_bending_bottom = (model.moment / EI) * (building_height / 2)

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
