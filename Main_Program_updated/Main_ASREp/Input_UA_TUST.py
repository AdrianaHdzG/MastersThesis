import sys
sys.path.append('ASREpy-main')
sys.path.append('FunctionScripts')

# %% LOAD TUNNEL PARAMETERS
z_t              = 8.00                  # [m] Tunnel center depth
vlt             = 0.01                  # [-] Relative Tunnel volume loss
Kt              = 0.5                   # [-] Estimated
D_t             = 6.00                  # [m] Diameter of tunnel

# %% LOAD BUILDING INFORMATION
building_height = 9                     # [m] Building height
building_width  = 0.5
Eb              = 1e9
poissons_ratio  = 0.3

neutral_line    = 4.5                   # [m] Neutral axis IN PURE bending. # Distance from BOTTOM FIBRE of BEAM TO NEUTRAL AXIS
dist_a          = 4.5                   # [m] Location of beam-soil interface. 0 = beam axis, h/2 beam bottom fibre

# %% Soil properties
Es              = 100                   # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)
soil_poisson    = 0.5                   # [-]


