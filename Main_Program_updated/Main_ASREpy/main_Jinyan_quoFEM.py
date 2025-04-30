# %% Imports
import numpy as np

from Input_general import *  # import * = import all
from Wall_deflection_v_38 import wall_deflection
from strain_analysis_greenfield import *
from FEA_LTS import categorize_damage, compute_tensile_strain_Jinyan, generate_output
from localStiff3D import *
from Output_uncertainty_analysis import *
import sys
from prepare_greenfield_disp import *
# from plotDisp import *

plot_var = "yes"  # "yes" or "no" to define if plot wanted or not


ASREpy_dir = os.path.join(os.path.dirname(__file__), "..", "ASREpy-main")
sys.path.append(ASREpy_dir)
sys.path.append("D:\Main_Program_updated\ASREpy-main")
sys.path.append("S:\Main_Program_updated\ASREpy-main")
import ASREpy

if avg_wall_disp_construction <= 0:
    avg_wall_disp_construction = 0
if avg_wall_disp_installation <= 0:
    avg_wall_disp_installation = 0
if avg_wall_disp_construction == 0 and avg_wall_disp_installation == 0:
    if output == 'quoFEM':
        generate_output(output_fields, 0, output_fields, error_flag=1)
        sys.exit()
    elif output == 'OLDquoFEM':
        with open("results.out", "w") as f:
            f.write("{} {} {} {} {} {} {}".format(0, 0, 0, 0, 0, 0, 1))
        sys.exit()
    else:
        print('Error: output field must be either "quoFEM" or "OLDquoFEM')
        sys.exit()



# %% INSTALLATION EFFECTS
average_wall_displacement_normalized = avg_wall_disp_installation  # [-] Average wall displacement normalized by the wall height (as pure number)
shape_wall_deflection = shape_wall_deflection_i  # [-] Shape of wall deflection M0-M4 | 0 = Uniform, 1 = Cantilever,
#     2 = Parabola type, 3 = Composite type, 4 = Kick-in type, 5 custom

if average_wall_displacement_normalized <= 0:  # Do not allow values below 0
    print(
        f"Average wall displacement (installation) is {average_wall_displacement_normalized} and is therefore set to 0")
    average_wall_displacement_normalized = 0
    vertical_displacement_installation = np.zeros(
        [1, num_nodes])  # Installation effects vertical settlement at building
    horizontal_displacement_installation = np.zeros(
        [1, num_nodes])  # Installation effects horizontal movement at building
else:
    Sx, Sz, Deltaw_all, VLW, Deltaw_max, check_ratio, VLS, ratiovol = wall_deflection(retaining_wall_depth,
                                                                                      average_wall_displacement_normalized,
                                                                                      shape_wall_deflection_i,
                                                                                      soil_ovalization, volumetric,
                                                                                      x_coords=building_coords,
                                                                                      C1_val=C1, C2_val=C2)

    vertical_displacement_installation = Sz  # Installation effects vertical settlement at building
    horizontal_displacement_installation = Sx  # Installation effects horizontal movement at building

# DISPLACEMENTS ARE COMBINED VIA SUPERPOSITION (ELASTIC SOIL BEHAVIOUR ASSUMED)
combined_vertical_displacement = vertical_displacement_installation
combined_horizontal_displacement = horizontal_displacement_installation

# If we want to include the plot, just call the wall_deflection function here with the full x span
"""
Sx, Sz, Deltaw_all, VLW, Deltaw_max, check_ratio, VLS, ratiovol = wall_deflection(retaining_wall_depth,
                                                                                  average_wall_displacement_normalized,
                                                                                  shape_wall_deflection,
                                                                                  soil_ovalization, volumetric,
                                                                                  C1, C2,
                                                                                  x_coords = input)
delta_v_x_ins_plot = Sx    # Installation effects 
delta_v_z_ins_plot = Sz    # Installation effects
"""

# %% Greenfield Settlements based on Wall Deflection Shape Functions
average_wall_displacement_normalized = avg_wall_disp_construction  # [-] Average wall displacement normalized by the wall height (as pure number)
shape_wall_deflection = shape_wall_deflection_c  # [-] Shape of wall deflection M0-M4 | 0 = Uniform, 1 = Cantilever,
#     2 = Parabola type, 3 = Composite type, 4 = Kick-in type, 5 custom

if average_wall_displacement_normalized <= 0:  # Do not allow values below 0
    print(
        f"Average wall displacement (construction) is {average_wall_displacement_normalized} and is therefore set to 0")
    average_wall_displacement_normalized = 0
    vertical_displacement_construction = np.zeros(
        [1, num_nodes])  # Installation effects vertical settlement at building
    horizontal_displacement_construction = np.zeros(
        [1, num_nodes])  # Installation effects horizontal movement at building
else:
    Sx, Sz, Deltaw_all, VLW, Deltaw_max, check_ratio, VLS, ratiovol = wall_deflection(retaining_wall_depth,
                                                                                      average_wall_displacement_normalized,
                                                                                      shape_wall_deflection,
                                                                                      soil_ovalization, volumetric,
                                                                                      x_coords=building_coords,
                                                                                      C1_val=C1, C2_val=C2)
    vertical_displacement_construction = Sz  # Installation effects vertical settlement at building
    horizontal_displacement_construction = Sx  # Installation effects horizontal movement at building

# DISPLACEMENTS ARE COMBINED VIA SUPERPOSITION (ELASTIC SOIL BEHAVIOUR ASSUMED)
combined_vertical_displacement += vertical_displacement_construction
combined_horizontal_displacement += horizontal_displacement_construction

# Rewrite here:
vertical_displacement_ground_building = combined_vertical_displacement
horizontal_displacement_ground_building = np.zeros([1, np.size(combined_horizontal_displacement, 1)])  # [m]
horizontal_displacement_ground_building = combined_horizontal_displacement

# %% LIMITING TENSILE STRAIN | GREENFIELD

# Green-field analysis for Limiting Tensile Strain method:
eps_tensile, max_angular_distortion_gf, max_eps_h_gf = strain_analysis_greenfield(vertical_displacement_ground_building,
                                                                                  horizontal_displacement_ground_building,
                                                                                  length_beam_element, length_beam,
                                                                                  poissons_ratio)

print("Tensile strain calculation for Greenfield analysis")
highest_damage_greenfield, max_tensile_eps_gf = categorize_damage(max(eps_tensile))  # Highest damage category

# %% Jinyan model displacements
# Horizontal and vertical displacements gone through interpolation
# combined_vertical_displacement      = np.array(combined_vertical_displacement).flatten() * 100
# combined_horizontal_displacement    = np.array(combined_horizontal_displacement).flatten() * 100
combined_vertical_displacement = -np.array(vertical_displacement_ground_building).flatten()  # STD: FEM CONVENTION
combined_horizontal_displacement = np.array(horizontal_displacement_ground_building).flatten()

# UPWARDS POSITIVE, TOWARDS RIGHT POSITIVE AND COUNTERCLOCKWISE POSITIVE


## Jinyan model inputs
beamY = np.zeros_like(building_coords)
beamZ = np.zeros_like(building_coords)
qfoot = 0  # [N/m] The uniform weight applied on the beam. Unit: N/m (Weight along unit length in the longitudinal direction)


## Jinyan / Peter & Seb modified elastic model
model_el = ASREpy.ASRE_Timoshenko_model_ps(building_coords.size, building_coords, beamY, beamZ, building_height,
                                        building_width, solver='elastic')
model_el.set_beam_properties(Eb, Eb / Gb, qfoot, d_a=dist_a)
model_el.set_soil_properties(Es_isotropic, soil_poisson, mu_int=0)
'''
# Currently use standard model:
model_el = ASREpy.ASRE_Timoshenko_model(building_coords.size, building_coords, beamY, beamZ, building_height,
                                        building_width, solver='elastic')
model_el.set_beam_properties(Eb, Eb / Gb, qfoot, d_NA=dist_a)
model_el.set_soil_properties(Es_isotropic, soil_poisson, mu_int=0)
'''

model_el.run_model(combined_horizontal_displacement, np.zeros_like(combined_horizontal_displacement),
                   combined_vertical_displacement, 'strain+disp+force')

model_properties = np.array([EA, EI, Gb, Ab, poissons_ratio, shear_factor_midpoint, building_height,
                             dist_a, neutral_line])

# %% Convert forces if d_NA != 0
# print('before conversion', model_el.axialForce)
if dist_a != 0:
    model_el, total_disp = moveResultLocation(model_el, dist_a, num_nodes)
    F_M_deltaT_el, F_N_deltaT_el, F_S_deltaT_el = compute_internal_forces(EA, EI, EI, GAs, GAs, length_beam_element,
                                                                          Gb, 1, total_disp, num_nodes)
    # Use updated axial force
    model_el.axialForce = F_N_deltaT_el


# print('after conversion', model_el.axialForce)
'''
# %% Jinyan model forces
data = {  # Make a DATA struct for all models run
    'model_elastic': {
        'x_coordinate': building_coords,
        'beam_strain_top': model_el.beam_strain_top,
        'beam_strain_bottom': model_el.beam_strain_bottom,
        'beam_strain_diagonal': model_el.beam_strain_diagonal,
        'beam_strain_axial': model_el.axialForce / EA,
        'moment': model_el.moment,
        'axialForce': model_el.axialForce,
        'shearForce': model_el.shearForce,
        'dispVertical': model_el.beam_DispV,  # z coordinate (Displacements at soil surface)
        'dispLongitudinal': model_el.beam_DispL  # x coordinate (Displacements at soil surface)
        # Note to self, recalculate displacements at beam axis.
    }
}
'''

# %% Recalculate strains to ensure strains are correct
dataReturn = compute_tensile_strain_Jinyan(model_el, model_properties)  # Call function to calculate new strains

# Calculate maximum tensile strains
max_eps_t_ssi = categorize_damage(dataReturn['max_e_t'])

# %% Validation done

# Keep strains without units [-]
if output == 'quoFEM':
    generate_output(dataReturn, max_eps_t_ssi, output_fields, error_flag=0)

elif output == 'OLDquoFEM':
    with open("results.out", "w") as f:
        f.write("{} {} {} {} {} {} {}".format(max(dataReturn['axial_strain_exx']), dataReturn['e_t_b'],
                                           dataReturn['max_e_xy'], dataReturn['max_e_t_mid'], dataReturn['max_e_t'],
                                           max_eps_t_ssi[0], 0))
else:
    print('this->')
