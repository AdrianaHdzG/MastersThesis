from Input_UA_TUST import *

# %% GENERIC INPUT
output = 'OLDquoFEM'   # 'quoFEM', 'OLDquoFEM' or 'General'
mode = 'SSI'           # 'SA' = Greenfield strain analysis, 'SSI' Soil-Structure-Interaction analysis. SSI includes SA
input_type = 'TUNNEL'                # 'TUNNEL' is the tunnelling case from Franza et al [2020], 'WALL' is default
solver = 'EL'                   # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver
# solver = 'EP'                   # 'EL' is cauchy elastic solver, 'EP' is elastoplasic solver

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