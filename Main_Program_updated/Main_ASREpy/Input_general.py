from Input_uncertainty_analysis import *
import numpy as np
# import scipy as sci

# %% INPUT PARAMETERS

output = 'quoFEM'  # 'quoFEM' or 'General'

# MATERIAL PARAMETERS
# Young_modulus       = Ebeam         # [Pa] Young's modulus
poissons_ratio      = 0.3           # [-]
coeff_tim           = 1             # 0 for Bernoulli-Euler and 1 for Timoshenko beam

# BUILDING PARAMETERS
building_offset     = 20             # [m] Offset from wall
building_width      = 1             # [m] Width of building foundation
length_beam         = 20            # [m] Length of beam (or building)
Ab                  = building_height * building_width  # [m^2] Cross-sectional building area
Eb                  = EA / Ab       # [Pa] Young's modulus of building
Ib                  = EI / Eb       # [m^4] Second moment of area 
Gb                  = GAs / (Ab * (10 + 10 * poissons_ratio) / (12 + 11 * poissons_ratio))  # [Pa] Shear modulus


# SOIL PARAMETERS
soil_poisson        = 0.3           # [-]
includeVesic        = False         # [True/False] Allow for vesic calibration of winkler springs
soil_ovalization    = 0             # [-] Soil ovalization
volumetric          = 1             # [-] Volumetric

# RETAINING WALL PARAMETERS
retaining_wall_depth        = 25    # [m] Depth of the retaining wall
shape_wall_deflection_i     = 2     # [-] Shape of wall deflection for installation effects 
shape_wall_deflection_c     = 5     # [-] Shape of wall deflection for construction/excavation effects 
                                    # Shape of wall deflection M0-M4 | 0 = Uniform, 1 = Cantilever, 2 = Parabola type, 3 = Composite type, 4 = Kick-in type, 5 Custom

# BUILDING DAMAGE
trough_width_param              = 0.5   # [-] Trough width parameter (approx 0.5 for fine-grained, 0.3 for coarse)
vertical_to_horizontal_gradient = 1.0   # [-] zv/zt - vertical to horizontal gradient


# %% MESH PARAMETERS

# GEOMETRY
num_nodes                 = 21     # Number of nodes
num_elements              = num_nodes - 1         # Number of elements
degrees_freedom_per_node  = 3                     # Degrees of freedom per node
num_global_DOF            = degrees_freedom_per_node * num_nodes  # Total number of DEGREES OF FREEDOM

# SUPPORT AND SOIL NODES
# support_nodes = np.array([0, num_nodes - 1])    # Supports
soil_nodes = np.arange(0, num_nodes)            # Winkler springs 

# BUILDING COORDINATES
building_coords           = building_offset + np.linspace(0, length_beam, num_nodes)  # x-coordinates for the building [m]

# %% SPRING CALIBRATION 

# SPRING STIFFNESS CALIBRATION
stiffness_vertical_spring    *= 1E6         # Rewrite to N/m
stiffness_horizontal_spring  *= 1E6
Es_isotropic                  = stiffness_vertical_spring  # Isotropic soil
rotational_stiffness_spring   = 0 * 1E6


# %% PARAMETER CALIBRATION AND CONVERSION

# WALL DISPLACEMENT CONVERSION TO PERCENT
avg_wall_disp_installation = avg_wall_disp_installation / 100    # Installation effects
avg_wall_disp_construction = avg_wall_disp_construction / 100       # Construction sequence

# CONVERSION OF PARAMETERS (AUTOMATIC RECTANGULAR CROSSSECTION)
total_degrees_freedom     = num_nodes * degrees_freedom_per_node         # Total degrees of freedom
length_beam_element       = length_beam / num_elements                   # Length of beam element [m]
shear_factor_midpoint     = 1.5                                          # Highest recorded shear stress at beam midpoint

# CONTROLLING THE WALL CALCULATION
coords_normal_wall        = np.arange(length_beam_element, 100 + length_beam_element, length_beam_element)  # [m] Coordinates in the direction normal to the wall
average_wall_displacement_normalized = avg_wall_disp_installation  # [-] Average wall displacement normalized by the wall height

# %% DRAWDOWN RELATED INPUT

# GROUNDWATER MOVEMENT / CORRECTION
well_diameter             = 0.5                 # Diameter of well
depth_of_well             = 30                  # Height of water table to bottom of aquifer [m]
drawdown                  = 0                   # Drawdown at the well [m]
# drawdown                  = drawdown_size       # Drawdown at the well [m]
influence_radius          = 250                 # Influence radius of pump [m]
specific_weight_water     = 9.82e3              # Specific weight of water [N/m^3]

# SOIL PARAMETERS
specific_weight_unsat_soil = 17e3                 # Specific weight of unsaturated soil [N/m^3]
specific_weight_sat_soil  = 20e3                  # Specific weight of saturated soil [N/m^3]
compression_index         = 0.3                   # Compression index [-]
initial_void_ratio        = initial_void_ratio_param       # Initial void ratio [-]
watertable_depth_z        = 2
height_soil_layer         = depth_of_well + watertable_depth_z     # Height of whole soil layer [m]


# UNIT CONVERSION, 1 for [N], 1000 for output [kN] etc.
# convert_units = 1000
