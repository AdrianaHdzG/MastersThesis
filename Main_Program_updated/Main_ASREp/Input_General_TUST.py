from Input_UA_TUST import *
import numpy as np


# %% LOAD TUNNEL PARAMETERS
z_t              = 8.00                  # [m] Tunnel center depth
# vlt             = 0.01                  # [-] Relative Tunnel volume loss
# Kt              = 0.5                   # [-] Estimated
D_t             = 6.00                  # [m] Diameter of tunnel

# %% LOAD BUILDING INFORMATION
building_height = 9                     # [m] Building height
building_width  = 0.5
# Eb              = 1e9
# poissons_ratio  = 0.3

neutral_line    = 4.5                   # [m] Neutral axis IN PURE bending. # Distance from BOTTOM FIBRE of BEAM TO NEUTRAL AXIS
dist_a          = 4.5                   # [m] Location of beam-soil interface. 0 = beam axis, h/2 beam bottom fibre

# %% Soil properties
# Es              = 100                   # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)
soil_poisson    = 0.5                   # [-]


# %% GENERIC INPUT
output = 'OLDquoFEM'   # 'quoFEM', 'OLDquoFEM' or 'General'
mode = 'SSI'           # 'SA' = Greenfield strain analysis, 'SSI' Soil-Structure-Interaction analysis. SSI includes SA
def_mode = 'Shear'     # Primary Deformation mode of greenfield structure, 'Shear' or 'Bending'
input_type = 'TUNNEL'  # 'TUNNEL' is the tunnelling case from Franza et al [2020], 'WALL' is default
solver = 'EL'          # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver
# solver = 'EP'        # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver

# %% GEOMETRY
length_beam = 20               # Beam length meters
num_nodes = 41                 # Number of nodes
num_elements = num_nodes - 1
eccentricity = 0.5             # 0 or 0.5


# %% BUILDING PARAMETERS
Ib = 1 / 12 * building_height**3 * building_width
Ab = building_height * building_width
As = Ab * (10 + 10 * poissons_ratio) / (12 + 11 * poissons_ratio)
Gb = Eb / (2 + 2 * poissons_ratio)
length_beam_element = length_beam / num_elements

EA = Eb * Ab
EI = Eb * Ib
GAs = Gb * As


# %% Update soil parameters
Es_isotropic                  = Es * 1E6  # [Pa] Isotropic soil

# %% GENERIC
shear_factor_midpoint = 1.5

# %% FOR ELASTOPLASTIC ANALYSIS
mu = 30  # [deg] Friction angle in degrees
mu_int = np.tan(mu * np.pi / 180)  # Calculation to rad


# LOADS ON FOOTING
qfoot = 3.2 * 10 * 1000  # [N] Load on the footing per running meter, necessary for good results
