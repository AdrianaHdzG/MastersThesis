# This script defines the uncertainty variables used in quoFEM

# MATERIAL PARAMETERS

EA      = 3.06E9                        # [N]
GAs     = 8.03E8                        # [N]
EI      = 2.30E9                        # [N*m^2]

# BUILDING INFORMATION
building_height = 3.00                   # [m] From building top to bottom foundation
neutral_line    = 1.00                   # [m] Neutral axis IN PURE bending. # Distance from BOTTOM FIBRE of BEAM
dist_a          = 0.00                   # [m] Location of beam-soil interface. 0 = beam axis, h/2 beam bottom fibre
# level as soil. Any value offset, where neutral axis is above ground level is treated as positive.


# WALL DISPLACEMENT
avg_wall_disp_installation  = 0.04      # [%] Installation effects (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.029244)
avg_wall_disp_construction  = 0.38      # [%] Construction sequence (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.7311)
C1                          = 0.30      # [-]
C2                          = 0.30      # [-]

# SOIL INITIAL PARAMETERS
Es                          = 27        # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)



