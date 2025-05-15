import sys
import numpy as np
sys.path.append('ASREpy-main')
sys.path.append('FunctionScripts')
import numpy as np
import os
import sys

# from plotFunctions import plot_disp, plot_strain
# from strain_cal import compute_tensile_strain_Jinyan
# from localStiff3D import *

# %% LOAD TUNNEL PARAMETERS
z0 = 8.00           # [m] Tunnel center depth
vlt = 1.0 / 100     # [-] Relative Tunnel volume loss
Kt = 0.5            # [-] Estimated
d = 6.00            # [m] Diameter of tunnel

# %% LOAD BUILDING INFORMATION
building_height = 9                                     # [m] Building height
building_width = 0.5
Eb = 1e9
poissons_ratio = 0.3

neutral_line    = building_height / 2                   # [m] Neutral axis IN PURE bending. # Distance from BOTTOM FIBRE of BEAM TO NEUTRAL AXIS
dist_a          = 0  # 4.5                                     # [m] Location of beam-soil interface. 0 = beam axis, h/2 beam bottom fibre

# %% Soil properties
Es   = 100                                              # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)
soil_poisson        = 0.5           # [-]


