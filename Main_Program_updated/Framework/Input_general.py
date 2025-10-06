from Input_uncertainty_analysis import *
import numpy as np
import scipy


# %% GENERAL INPUT PARAMETERS, NEEDED FOR GREENFIELD AND SSI


output = 'OLDquoFEM'   # 'quoFEM', 'OLDquoFEM' or 'General'
mode = 'SSI'           # 'SA' = Greenfield strain analysis, 'SSI' Soil-Structure-Interaction analysis. SSI includes SA
def_mode = 'Shear'     # Primary Deformation mode of greenfield structure, 'Shear' or 'Bending'


integration_mode = 'CValues'       # 'CValues' uses C1, C2, C3 or modes. 'Direct' uses a datafile like below.
# integration_mode = 'Direct'
input_type = 'WALL'                # 'TUNNEL' is the tunnelling case from Franza et al [2020], 'WALL' is default. "3DWALL" to use the 3D greenfield
solver = 'EL'                      # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver

if integration_mode == 'Direct':  # Example of how to load data
    data_wall = scipy.io.loadmat("FunctionScripts\example_wall_data.mat")


# BUILDING PARAMETERS
building_height     = 36.0          # [m] From building top to bottom foundation
neutral_line        = 18.0          # [m] Neutral axis IN PURE bending. # Distance from BOTTOM FIBRE of BEAM
building_offset     = 11.0          # [m] Offset from wall
building_width      = 1             # [m] Width of building foundation
length_beam         = 12            # [m] Length of beam (or building)
foundation_depth    = 2.2           # [m] Depth positive 2.2


# SOIL PARAMETERS
soil_poisson        = 0.5           # [-] Soil poisson ratio, recommended 0.5 for undrained and 0.2-0.3 for drained
soil_ovalization    = 0             # [-] Soil ovalization
volumetric          = 1             # [-] Volumetric

# RETAINING WALL PARAMETERS
# retaining_wall_depth        = 20.7  # [m] Depth of the retaining wall
C1                          = -0.36  # [-] day 114 -0.36
C2                          = 1.29    # [-] 1.29 
shape_wall_deflection_i     = 2       # [-] Shape of wall deflection for installation effects
shape_wall_deflection_c     = 5       # [-] Shape of wall deflection for construction/excavation effects



# % MESH PARAMETERS
# GEOMETRY
num_nodes                 = 101     # Number of nodes
num_elements              = num_nodes - 1         # Number of elements

# BUILDING COORDINATES
building_coords           = building_offset + np.linspace(0, length_beam, num_nodes)  # x-coordinates for the building [m]


# % PARAMETER CALIBRATION AND CONVERSION

# WALL DISPLACEMENT CONVERSION TO PERCENT
avg_wall_disp_installation = 0  # avg_wall_disp_installation / 100    # Installation effects
avg_wall_disp_construction = avg_wall_disp / 100       # Construction sequence

# Calculate element length
length_beam_element       = length_beam / num_elements                   # Length of beam element [m]

# % 3D GREENFIELD PARAMETERS
# These parameters are only used if integration_mode = '3DWall'
# Station box dimensions
L_x_3D = 10          # [m] Half-width of the station box
L_y_3D = 30          # [m] Half-length of the station box
He_Hwratio_3D = 1.0               # [-] Ratio of He/Hw

y0_line         = 0.0    # [m] fixed y for the line at y = 0 centerline of the station box representative of 2D case
z0_line         = foundation_depth  # [m] elevation of the line (0=ground surface)

# Avg_wall_displacement (beta) for each wall
beta_CCS_wall_1_3D = avg_wall_disp / 100
beta_CCS_wall_2_3D = avg_wall_disp / 100
beta_CCS_wall_3_3D = avg_wall_disp / 100
beta_CCS_wall_4_3D = avg_wall_disp / 100

# DMM parameters for 3D wall deflection shape (match MATLAB C1, C2 values)
C1_3D                       = C1    # [-] day 114
C2_3D                       = C2    # [-]

building_offset_3D = building_offset  # [m] Offset from wall
building_coords_3D = (L_x_3D + building_offset) + np.linspace(0.0, length_beam, num_nodes)
length_beam_3D = length_beam  # [m] Length of beam (or building)

# Vertical wall deflection shape (match MATLAB switches)
#   3  : Parabolic
#   30 : Parabolic + Mu & Huang longitudinal reduction (Gaussian)
#   31 : Parabolic + Roboski & Finno (erfc-type)
#   5  : DMM combo (C1 cantilever, C2 parabolic, C3 kick-in)
#   50 : DMM + Mu & Huang
#   51 : DMM + Roboski & Finno
switch_shape_3D = 5

# Discretization for the cavity stacks
delta_z_cavities_3D   = retaining_wall_depth / retaining_wall_depth   # vertical bins (m)
delta_xyperimeter_3D  = 2.5                           # along-wall segment size (m)

# Solution type for symmetry/taper:
#   1 = single wall mirrored
#   2 = 4 walls analytical (no taper)
#   3 = 4 walls semi-analytical (with taper)  ← recommended
switch_solution_type_3D = 3



# %% IF MODE = SSI
if mode == 'SSI':
    # MATERIAL PARAMETERS Conversion

    EA = EA * 1E9  # [N]
    GAs = GAs * 1E9  # [N]
    EI = EI * 1E9  # [N*m^2]

    Ab = 2.76  # [m^2] Cross-sectional building area
    Eb = EA / Ab  # [Pa] Young's modulus of building
    Ib = EI / Eb  # [m^4] Second moment of area
    Gb = GAs / 1.65  # [Pa] Shear modulus

    # BUILDING INFORMATION
    dist_a = 18.0  # [m] Location of beam-soil interface. 0 = beam axis, h/2 beam bottom fibre
    # level as soil. Any value offset, where neutral axis is above ground level is treated as positive.

    # % SPRING CALIBRATION
    # SPRING STIFFNESS CALIBRATION
    Es_isotropic = Es * 1E6  # Isotropic soil (N/m/m)

    # LOADS ON FOOTING
    qfoot = 3.2 * 10 * 1000  # [N] Load on the footing per running meter, necessary for good results

    # Shear factor
    shear_factor_midpoint = 1.5  # Highest recorded shear stress at beam midpoint


    if solver == 'EP':
        mu_int = np.tan(30 * np.pi / 180)



