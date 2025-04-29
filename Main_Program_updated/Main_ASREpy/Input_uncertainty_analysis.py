# This script defines the uncertainty variables used in quoFEM

# MATERIAL PARAMETERS
EA      = 3.06E9                        # [N]
GAs     = 8.03E8                        # [N]
EI      = 2.30E9                        # [N*m^2]

# BUILDING INFORMATION
# neutral_line    = 0.01                  # [m] Distance from bottom of foundation to neutral-line
building_height = 3                     # [m] From building top to bottom foundation
dist_NA         = 0                     # [m] Like neutral-line but for Jinyan analysis. 0 indicates dist_NA is same
# level as soil. Any value offset, where neutral axis is above ground level is treated as positive.


# WALL DISPLACEMENT
avg_wall_disp_installation  = 0.04      # [%] Installation effects (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.029244)
avg_wall_disp_construction  = 0.38    # [%] Construction sequence (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.7311)
C1                          = 0.30     # [-]
C2                          = 0.30      # [-]
# avg_wall_disp_installation  = -0.04      # [%] Installation effects (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.029244)
# avg_wall_disp_construction  = -0.38    # [%] Construction sequence (COV = 73.11% according to Zhao and DeJong, 2023 --> standard deviation = 0.7311)
# SOIL INITIAL PARAMETERS
initial_void_ratio_param    = 1         # [-] Initial_void_ratio being 1.0
stiffness_vertical_spring   = 27        # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)
stiffness_horizontal_spring = 27        # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023.

# GROUNDWATER MOVEMENT / CORRECTION
# drawdown_size   = 0.0                   # [m] Drawdown at the well

