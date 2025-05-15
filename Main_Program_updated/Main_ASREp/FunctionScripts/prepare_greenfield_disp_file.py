from LoadPackages import *


def prepare_greenfield_disp_file(val_data_dir, vl, beam_id, beamX):
    ## Prepare Ux-GF
    # Read Ux-GF and Uz-GF
    dispX_gf_file = os.path.join(val_data_dir, f'FranzaDeJong{beam_id}-'
                                               f"VL{str(vl * 100).replace('.', 'p')}-Ux-GF.csv")
    dispX_gf = pd.read_csv(dispX_gf_file, header=None)

    dispZ_gf_file = os.path.join(val_data_dir, f'FranzaDeJong{beam_id}-'
                                               f"VL{str(vl * 100).replace('.', 'p')}-Uz-GF.csv")
    dispZ_gf = pd.read_csv(dispZ_gf_file, header=None)
    # Add the left half of Ux-GF and Uz-GF
    left_half = dispX_gf.copy()
    left_half[0] = -left_half[0]
    left_half[1] = -left_half[1]
    dispX_gf = pd.concat([left_half, dispX_gf])
    dispX_gf = dispX_gf.sort_values(by=0).reset_index(drop=True)

    left_half = dispZ_gf.copy()
    left_half[0] = -left_half[0]
    dispZ_gf = pd.concat([left_half, dispZ_gf])
    dispZ_gf = dispZ_gf.sort_values(by=0).reset_index(drop=True)
    # Interpolate the dispX_gf with beamX mesh
    f = interpolate.interp1d(dispX_gf[0], dispX_gf[1], fill_value='extrapolate')
    dispX_gf = f(beamX)

    f = interpolate.interp1d(dispZ_gf[0], dispZ_gf[1], fill_value='extrapolate')
    dispZ_gf = f(beamX)
    return dispX_gf / 1000, dispZ_gf / 1000


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


def plotSimple(x_coord, yval, label):
    fig, ax = plt.subplots()
    ax.plot(x_coord, yval * 1000)
    ax.set_xlabel('x_coord')
    ax.set_ylabel(label)
    ax.set_title(f'Plot of X-coordinate [m] vs {label} [mm]')
    plt.show()



