# This script defines the uncertainty variables used in quoFEM

# MATERIAL PARAMETERS

EA      = 25.9                        # [GN]
GAs     = 5.95                        # [GN]
EI      = 2777                        # [GN*m^2]

# MATERIAL PARAMETERS
poissons_ratio      = 0.3           # [-] For the building materials

# Wall parameters
avg_wall_disp               = 0.33      # [%] day 114, Construction sequence (COV = 73.11% for tunnelling according to Zhao and DeJong, 2023) Same COV considered here
retaining_wall_depth        = 20.7      # [m] Depth of the retaining wall from surface to rigid layer

# SOIL INITIAL PARAMETERS

Es                          = 175        # [MN/m] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14)


