from Input_uncertainty_analysis import *
import numpy as np
import scipy

# %% INPUT PARAMETERS


output = 'OLDquoFEM'   # 'quoFEM', 'OLDquoFEM' or 'General'
mode = 'SSI'           # 'SA' = Greenfield strain analysis, 'SSI' Soil-Structure-Interaction analysis. SSI includes SA
integration_mode = 'CValues'       # 'CValues' uses C1, C2, C3 or modes. 'Direct' uses a datafile like below.
# integration_mode = 'Direct'
input_type = 'WALL'                # 'TUNNEL' is the tunnelling case from Franza et al [2020], 'WALL' is default
solver = 'EL'                      # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver

if integration_mode == 'Direct':  # Example of how to load data
    data_wall = scipy.io.loadmat("FunctionScripts\example_wall_data.mat")

# MATERIAL PARAMETERS
poissons_ratio      = 0.3           # [-] For the building materials

# BUILDING PARAMETERS
building_offset     = 20            # [m] Offset from wall
building_width      = 1             # [m] Width of building foundation
length_beam         = 20            # [m] Length of beam (or building)
Ab                  = building_height * building_width  # [m^2] Cross-sectional building area
Eb                  = EA / Ab       # [Pa] Young's modulus of building
Ib                  = EI / Eb       # [m^4] Second moment of area
Gb                  = GAs / (Ab * (10 + 10 * poissons_ratio) / (12 + 11 * poissons_ratio))  # [Pa] Shear modulus
foundation_depth    = 2.2           # [m] Depth positive


# SOIL PARAMETERS
soil_poisson        = 0.5           # [-] Soil poisson ratio, recommended 0.5 for undrained and 0.2-0.3 for drained
soil_ovalization    = 0             # [-] Soil ovalization
volumetric          = 1             # [-] Volumetric

# RETAINING WALL PARAMETERS
retaining_wall_depth        = 25.0  # [m] Depth of the retaining wall
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

# BUILDING COORDINATES
building_coords           = building_offset + np.linspace(0, length_beam, num_nodes)  # x-coordinates for the building [m]

# %% SPRING CALIBRATION
# SPRING STIFFNESS CALIBRATION
Es_isotropic                  = Es * 1E6  # Isotropic soil (N/m/m)

# LOADS ON FOOTING
qfoot = 3.2 * 10 * 1000  # [N] Load on the footing per running meter, necessary for good results
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

# %% FOR ELASTOPLASTIC ANALYSIS
if solver == 'EP':
    mu_int = np.tan(30 * np.pi / 180)



