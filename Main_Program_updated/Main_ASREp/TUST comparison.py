import sys

import numpy as np
sys.path.append('ASREpy-main')
sys.path.append('FunctionScripts')
sys.path.append('S:\Main_Program_updated\ASREpy-main')
from LoadPackages import *

from plotFunctions import plot_disp, plot_strain
from strain_cal import compute_tensile_strain_Jinyan
from localStiff3D import *
from prepare_greenfield_disp_file import prepare_greenfield

# %% Tunnel parameters
# beam_id = "STR-1"
z0 = 8  # Tunnel center depth[m]
vlt = 1.0 / 100  # Relative Tunnel volume loss
Kt = 0.5  # Estimated
d = 6  # Diameter of tunnel
l_b = 20  # Beam length meters
numNodes = 41  # Number of nodes
num_elements = numNodes - 1
H_b = 9  # Beam height
eccentricity = 0
d_na = H_b / 2  # Distance from normal line in beam to force application

# dispX_gf, dispZ_gf = prepare_greenfield_disp(val_data_dir, vlt, beam_id, footing_coord_x)

Uz, Ux, x_coord = prepare_greenfield(z0, vlt, Kt, d, l_b, eccentricity, numNodes)
# plotSimple(x_coord, Uz, 'Uz')
# plotSimple(x_coord, Ux, 'Ux')

plotfig7 = 'false'
plotfig8 = 'false'


# %% Create a frame model that is equivalent to STR-1
beamX = x_coord  # Using Jinyan notation, x coordinates of the problem

dfoot = H_b  # Height of footing or modelling
bfoot = 0.5
Ib = 1 / 12 * dfoot**3 * bfoot
As = dfoot * bfoot * (10 + 10 * 0.3) / (12 + 11 * 0.3)
beamY = np.zeros_like(beamX)
beamZ = np.zeros_like(beamX)
Eb = 1e9
EoverG = 2.6  # Assume very large G to model the Euler–Bernoulli beam used in Franza and DeJong
# EoverG = 0.00001

qfoot = 3.2 * 10 * 1000
Es = 100e6
nis = 0.5
# print(1/12 * dfoot**3 * bfoot * Eb / 1E9, 'EI') # Validate everything is implemented correctly
# print(dfoot * bfoot * Eb / 1E9, 'EA')
model_properties = np.array([Eb * dfoot * bfoot,  Eb * Ib, Eb / EoverG, dfoot * bfoot, 0.3, 1.5, dfoot, d_na])


    # %% Run elastoplastic model
mu_int = np.tan(30 * np.pi / 180)
model = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                         beamZ, dfoot, bfoot)
model.set_beam_properties(Eb, EoverG, qfoot, d_NA=0)
model.set_soil_properties(Es, nis, mu_int)

model.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')  # 'disp','strain+disp', 'strain+disp+force'

mu_int = np.tan(30 * np.pi / 180)
modelNA = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                           beamZ, dfoot, bfoot)
modelNA.set_beam_properties(Eb, EoverG, qfoot, d_NA=d_na)
modelNA.set_soil_properties(Es, nis, mu_int)

modelNA.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')  # 'disp','strain+disp', 'strain+disp+force'

# Correction to account for displacements being given at bottom node.
modelNA, total_disp = moveResultLocation(modelNA, d_na, numNodes)
F_M_deltaT_el1, F_N_deltaT_el1, F_S_deltaT_el1 = compute_internal_forces(Eb * dfoot * bfoot, Eb * Ib, Eb * Ib,
                                                                      Eb / EoverG * As, Eb / EoverG * As,
                                                                      l_b / (numNodes - 1), Eb / EoverG, 1, total_disp,
                                                                      numNodes)

modelNA.axialForce = F_N_deltaT_el1

data1 = compute_tensile_strain_Jinyan(modelNA, model_properties)


# %% Run elastic model

model_el = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                            beamZ, dfoot, bfoot,
                                            solver='elastic')
model_el.set_beam_properties(Eb, EoverG, qfoot, d_NA=0)
model_el.set_soil_properties(Es, nis, mu_int)

model_el.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')

# Validate data, as bending strains are wrongfully calculated atm.
data1 = compute_tensile_strain_Jinyan(model_el, model_properties)

model_elNA = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                              beamZ, dfoot, bfoot,
                                              solver='elastic')
model_elNA.set_beam_properties(Eb, EoverG, qfoot, d_NA=d_na)
model_elNA.set_soil_properties(Es, nis, mu_int)

model_elNA.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')


# %% Translate from outer fibre displacement into beam axis displacements
model_elNA, total_disp = moveResultLocation(model_elNA, d_na, numNodes)
F_M_deltaT_el, F_N_deltaT_el, F_S_deltaT_el = compute_internal_forces(Eb * dfoot * bfoot, Eb * Ib, Eb * Ib,
                                                                      Eb / EoverG * As, Eb / EoverG * As,
                                                                      l_b / (numNodes - 1), Eb / EoverG, 1, total_disp,
                                                                      numNodes)

model_elNA.axialForce = F_N_deltaT_el

data = compute_tensile_strain_Jinyan(model_elNA, model_properties)

print(data['exx,t,b'])

"""
internal_forces_vector = compute_internal_forces(Eb * dfoot * bfoot, Eb * Ib, Eb * Ib, Eb / EoverG * As,
                                                 Eb / EoverG * As, l_b / num_elements, Eb / EoverG, 1, total_disp,
                                                 num_elements)
axialForce = internal_forces_vector[:, 0]
shearForce = internal_forces_vector[:, 2]
moment = internal_forces_vector[:, 4]
"""
'''
# print(coorFactor)
# print(model_elNA.beam_DispL)
# print(model_elNA.beam_DispL)
# N_axis         = -1 * theta            * d_na

# rotationModified = np.zeros_like(model_elNA.moment)
# print('before conversion', model_elNA.axialForce)
# for i in range(model_elNA.moment.shape[0]):
#     rotationModified[i] = model_elNA.beam_RotaT[int((i + 1) / 2)]

# print(-1 * model_elNA.moment  / d_na)
# model_elNA.axialForce += -1 * model_elNA.moment / (Eb * 1 / 12 * dfoot**3 * bfoot) * d_na

# model_elNA.axialForce +=  (model_elNA.moment)  / (d_na * 1 / 3)
# print('After conversion', model_elNA.axialForce)
 #print('moment', model_elNA.moment)

from localStiff3D import*
#
degrees_freedom_per_node2 = 6
start_indices2 = np.arange(num_elements) * degrees_freedom_per_node2
end_indices2 = start_indices2 + 2 * degrees_freedom_per_node2
#
# Initialize element force matrix
element_force_matrix2 = np.zeros((2 * degrees_freedom_per_node2, num_elements))
print('element_force_matrix2.shape: ', element_force_matrix2.shape)
#
local_stiffness_matrix2 = localTim3D(Eb * dfoot * bfoot, Eb*Ib, Eb*Ib, Eb / EoverG * As, Eb / EoverG * As, l_b / num_elements, Eb / EoverG, 1) # localTim3D(EA, EIy, EIz, GAy, GAz, l, G, K): # ------------------------------------------------------------------------------------------------------- I added *100 --------------
#
total_disp = np.zeros([numNodes*degrees_freedom_per_node2,1])
#
# # u, v, w, delta_x, delta_y, delta_z
total_disp[0::6] = model_el.beam_DispL.reshape(-1, 1)      # [N] Horizontal forces
total_disp[1::6] = model_el.beam_DispT.reshape(-1, 1)      # [N] Translational Forces
total_disp[2::6] = model_el.beam_DispV.reshape(-1, 1)      # [N] Vertical forces
total_disp[3::6] = model_el.beam_RotaL.reshape(-1, 1)      # [Nm] Horizontal rotations
total_disp[4::6] = model_el.beam_RotaT.reshape(-1, 1)      # [Nm] Translational rotations
total_disp[5::6] = model_el.beam_RotaV.reshape(-1, 1)      # [Nm] Vertical rotations
#
# #print('stiff:' , local_stiffness_matrix2)
# #print('disp: ' , total_disp)
#
#
for start_idx, end_idx in zip(start_indices2, end_indices2):
    elem_force = (local_stiffness_matrix2 @ total_disp[start_idx:end_idx])
    element_force_matrix2[:, start_idx // degrees_freedom_per_node2] = elem_force.reshape(2 * degrees_freedom_per_node2)
#
print('element_force_matrix2.shape:', element_force_matrix2.shape)
#
element_force_matrix2_reduced = element_force_matrix2[[0, 2, 4, 6, 8, 10], :]
#
print('element_force_matrix2_reduced.shape', element_force_matrix2_reduced.shape)

'''




if plotfig7 == 'true':
    # %% Validaton plot
    plot_disp(beamX, Ux, 'Ux', beamX, model.beam_DispL, model_el.beam_DispL, modelNA.beam_DispL, model_elNA.beam_DispL, '7b')
    plot_disp(beamX, Uz, 'Uz', beamX, model.beam_DispV, model_el.beam_DispV, modelNA.beam_DispV, model_elNA.beam_DispV, '7a')


'''
# N_axis         = -1 * theta            * d_na
rotationModified = np.zeros_like(model_elNA.axialForce)
print('before conversion', model_elNA.axialForce)
for i in range(model_elNA.axialForce.shape[0]):
    rotationModified[i] = model_elNA.beam_RotaT[int((i + 1) / 2)]

#  M = E * I * theta
I_b = 1 / 12 * bfoot * dfoot**3
moment = Eb * I_b * rotationModified

# F = M / z         M = F * z
F = moment / d_na

#          N_0     = -F     + N_dna
model_elNA.axialForce = -1 * F + model_elNA.axialForce
print('after conversion', model_elNA.axialForce)
'''

# %% Save results for e/b = 0

import scipy.io

# Combine all data into a single dictionary
data = {
    'model_ep': {
        'x_coordinate': beamX,
        'beam_strain_top': model.beam_strain_top,
        'beam_strain_bottom': model.beam_strain_bottom,
        'beam_strain_diagonal': model.beam_strain_diagonal,
        'moment': model.moment,
        'axialForce': model.axialForce,
        'shearForce': model.shearForce,
    },
    'model_el': {
        'x_coordinate': beamX,
        'beam_strain_top': model_el.beam_strain_top,
        'beam_strain_bottom': model_el.beam_strain_bottom,
        'beam_strain_diagonal': model_el.beam_strain_diagonal,
        'moment': model_el.moment,
        'axialForce': model_el.axialForce,
        'shearForce': model_el.shearForce,
        'beam_strain_top1': data1['exx,t,b'],
        'beam_strain_bottom1': data1['exx,b,b'],
        'beam_strain_diagonal1': data1['tensile_strain_midpoint'],
    },
    'model_epNA': {
        'x_coordinate': beamX,
        'beam_strain_top': modelNA.beam_strain_top,
        'beam_strain_bottom': modelNA.beam_strain_bottom,
        'beam_strain_diagonal': modelNA.beam_strain_diagonal,
        'moment': modelNA.moment,
        'axialForce': modelNA.axialForce,
        'shearForce': modelNA.shearForce,
    },
    'model_elNA': {
        'x_coordinate': beamX,
        'beam_strain_top': model_elNA.beam_strain_top,
        'beam_strain_bottom': model_elNA.beam_strain_bottom,
        'beam_strain_diagonal': model_elNA.beam_strain_diagonal,
        'moment': model_elNA.moment,
        'axialForce': model_elNA.axialForce,
        'shearForce': model_elNA.shearForce,
    },
    'model_elNArecal': {
        'x_coordinate': beamX,
        'moment': model_elNA.moment,
        'axialForce': model_elNA.axialForce,
        'shearForce': model_elNA.shearForce,
        'beam_strain_top': data['exx,t,b'],
        'beam_strain_bottom': data['exx,b,b'],
        'beam_strain_diagonal': data['tensile_strain_midpoint'],
        'true_shear_strain': data['true_shear_strain'],
    },
    'model_epNARECAL': {
        'x_coordinate': beamX,
        'axialForce': F_N_deltaT_el1,
        'moment': modelNA.moment,
        'shearForce': modelNA.shearForce,
        'beam_strain_top': data1['exx,t,b'],
        'beam_strain_bottom': data1['exx,b,b'],
        'beam_strain_diagonal': data1['tensile_strain_midpoint'],
        'true_shear_strain': data1['true_shear_strain'],
    }
}

# print("data['model_elNArecal']['moment']", data['model_elNArecal']['moment'])

# Save the data to a .mat file
scipy.io.savemat("ELfig9NApy.mat", data)





    # %% Try new plotting rouitine
eccentricity = 0.5
d_na = H_b / 2  # Distance from normal line in beam to force application
Uz, Ux, x_coord = prepare_greenfield(z0, vlt, Kt, d, l_b, eccentricity, numNodes)
beamX = x_coord

    # %% Run elastoplastic model
mu_int = np.tan(30 * np.pi / 180)
model = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                         beamZ, dfoot, bfoot)
model.set_beam_properties(Eb, EoverG, qfoot, d_NA=0)
model.set_soil_properties(Es, nis, mu_int)

model.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')  # 'disp','strain+disp', 'strain+disp+force'

mu_int = np.tan(30 * np.pi / 180)
modelNA = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                           beamZ, dfoot, bfoot)
modelNA.set_beam_properties(Eb, EoverG, qfoot, d_NA=d_na)
modelNA.set_soil_properties(Es, nis, mu_int)

modelNA.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')  # 'disp','strain+disp', 'strain+disp+force'

# Correction to account for displacements being given at bottom node.
coorFactor = modelNA.beam_RotaT * d_na
modelNA.beam_DispL = coorFactor + modelNA.beam_DispL

    # %% Run elastic model

model_el = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                            beamZ, dfoot, bfoot,
                                            solver='elastic')
model_el.set_beam_properties(Eb, EoverG, qfoot, d_NA=0)
model_el.set_soil_properties(Es, nis, mu_int)

model_el.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')

model_elNA = ASREpy.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                              beamZ, dfoot, bfoot,
                                              solver='elastic')
model_elNA.set_beam_properties(Eb, EoverG, qfoot, d_NA=d_na)
model_elNA.set_soil_properties(Es, nis, mu_int)

model_elNA.run_model(Ux, np.zeros_like(Ux), Uz, 'strain+disp+force')

# Correction to account for displacements being given at bottom node.
coorFactor = model_elNA.beam_RotaT * d_na
model_elNA.beam_DispL = coorFactor + model_elNA.beam_DispL

if plotfig7 == 'true':
    # %% Validaton plot

    # plot_disp(beamX, Uz, 'Uz', beamX, model.beam_DispV, model_el.beam_DispV, modelNA.beam_DispV, model_elNA.beam_DispV, '7a')
    # plot_disp(beamX, Ux, 'Ux', beamX, model.beam_DispL, model_el.beam_DispL, modelNA.beam_DispL, model_elNA.beam_DispL, '7b')
    plot_disp(beamX, Uz, 'Uz', beamX, model.beam_DispV, model_el.beam_DispV, modelNA.beam_DispV, model_elNA.beam_DispV, '7c')
    plot_disp(beamX, Ux, 'Ux', beamX, model.beam_DispL, model_el.beam_DispL, modelNA.beam_DispL, model_elNA.beam_DispL, '7d')
    print('')



if plotfig8 == 'true':
    # F_N_deltaT_el_M, F_S_deltaT_el_M
    '''
    self.beam_strain_top = self.result_array_ptr[0:(self.nnode - 1) * 2]
    self.beam_strain_bottom = self.result_array_ptr[(self.nnode - 1) * 2:(self.nnode - 1) * 4]
    self.beam_strain_diagonal = self.result_array_ptr[(self.nnode - 1) * 4:(self.nnode - 1) * 6]
    result_array = self.result_array_ptr[(self.nnode - 1) * 6:(self.nnode - 1) * 6 + self.nnode * 6]
    self.beam_DispL = result_array[0::6]
    self.beam_DispT = result_array[1::6]
    self.beam_DispV = result_array[2::6]
    self.beam_RotaL = result_array[3::6]
    self.beam_RotaT = result_array[4::6]
    self.beam_RotaV = result_array[5::6]
    result_array = self.result_array_ptr[(self.nnode - 1) * 6 + self.nnode * 6:]
    self.moment = result_array[0:(self.nnode - 1) * 2]
    self.axialForce = result_array[(self.nnode - 1) * 2:(self.nnode - 1) * 4]
    self.shearForce = result_array[(self.nnode - 1) * 4:(self.nnode - 1) * 6]
    '''


    # print('true')
    # plot_strain(model_el, beamX)
print('done')



