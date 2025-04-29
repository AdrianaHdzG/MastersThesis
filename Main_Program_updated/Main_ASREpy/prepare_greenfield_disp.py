import numpy as np
import os
import pandas as pd
from scipy import interpolate

# %% GREENFIELD DISPLACEMENT | VECTOR FORMAT
def prepare_greenfield_disp_vec(beamX, disp_vec, direction):
    ## Prepare Ux-GF and Uz-GF with vectors
    
    # Stack vectors and convert to pandas DataFrame
    disp_gf = np.column_stack((beamX, disp_vec))
    disp_gf = pd.DataFrame(disp_gf)
    
    # Mirror the displacements about the z-axis
    left_half       = disp_gf.copy()
    left_half[0]    = -left_half[0]
    if direction == "horizontal":
        left_half[1]    = -left_half[1]     # Flip the horizontal displacement along the z-axis

    # Concatenate strings
    disp_gf = pd.concat([left_half, disp_gf])
    disp_gf = disp_gf.sort_values(by = 0).reset_index(drop = True)

    # Interpolate the disp_gf with beamX mesh
    f = interpolate.interp1d(disp_gf[0], disp_gf[1], fill_value='extrapolate')
    disp_gf = f(beamX)

    return disp_gf/1000


# %% GREENFIELD DISPLACEMENT | CSV FORMAT
def prepare_greenfield_disp_csv(beamX, csv_dir, filename, direction):
    ## Prepare greenfield displacements from csv file
    # beamX is the x-coordinates of the beam that we want interpolated
    # csv_dir is the directioary of the csv file
    # filename is the filename of the csv file (without ".csv") 

    # Find path and read csv
    disp_gf_file   = os.path.join(csv_dir, f"{filename}.csv")   
    disp_gf        = pd.read_csv(disp_gf_file, header=None)     

    # Mirror the displacements about the z-axis
    left_half       = disp_gf.copy()
    left_half[0]    = -left_half[0]
    if direction == "horizontal":
        left_half[1]    = -left_half[1]     # Flip the horizontal displacement along the z-axis

    # Concatenate strings
    disp_gf = pd.concat([left_half, disp_gf])
    disp_gf = disp_gf.sort_values(by = 0).reset_index(drop = True)

    # Interpolate the disp_gf with beamX mesh
    f = interpolate.interp1d(disp_gf[0], disp_gf[1], fill_value='extrapolate')
    disp_gf = f(beamX)

    return disp_gf/1000

def prepare_greenfield(z0, vlt, Kt, d, l_b, eccentricity, num):

    if eccentricity == 0:
        domainStart = -0.5 * l_b
        domainEnd = 0.5 * l_b
    if eccentricity == 0.5:
        domainStart = -l_b
        domainEnd = 0

    x_coord = np.linspace(domainStart, domainEnd, num)  # Create domain vec.

    i = Kt * z0  # Calculate horizontal distance to inflection point to tunnel centre line

    Uz_max = 1.25 * (d / 2) ** 2 * vlt / i  # Calculate maximum settlement

    Uz = np.zeros([num])
    Ux = np.zeros([num])
    for j in range(num):
        Uz[j] = -Uz_max * np.exp(-(x_coord[j] ** 2) / (2 * i ** 2))  # Calculate current settlement
        Ux[j] = Uz[j] * x_coord[j] / z0  # Calculate current deflection

    return Uz, Ux, x_coord