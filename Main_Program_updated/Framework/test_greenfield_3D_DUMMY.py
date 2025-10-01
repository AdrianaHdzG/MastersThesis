import sys, numpy as np
import matplotlib.pyplot as plt

sys.path.append('FunctionScripts')

from greenfield_3D_DMM import run_greenfield_3D_line
from strain_analysis import *

# --- Test inputs (match your typical station line to the right) ---
Hw = 20.0
L_x, L_y = 9.5, 32.0
building_offset = 11.0
length_beam = 12.0
num_nodes = 101
y0 =  0.0
z0 = -10.0

# --- Run 3D (tapered, SPy) ---
ux, uy, uz, zc, dw = run_greenfield_3D_line(
    Hw=Hw, L_x=L_x, L_y=L_y,
    switch_shape=5, C1=0.0, C2=1.0,           # parabolic only
    switch_solution_type=3,                   # semi-analytical taper
    building_offset=building_offset,
    length_beam=length_beam,
    num_nodes=num_nodes,
    y0=y0, z0=z0
)

print("Shapes (ux, uy, uz):", ux.shape, uy.shape, uz.shape)
print("Depth vec len:", len(zc), " Delta-wall vec len:", len(dw))
print("Max |ux|, |uy|, |uz|:", float(np.nanmax(np.abs(ux))), float(np.nanmax(np.abs(uy))), float(np.nanmax(np.abs(uz))))

# --- Rebuild x for plotting (same formula used internally) ---
x = (L_x + building_offset) + np.linspace(0.0, length_beam, num_nodes)

# --- Plots ---
plt.figure()
plt.plot(x, ux.ravel())
plt.title("Ux along building line")
plt.xlabel("x [m]"); plt.ylabel("Ux [m]"); plt.grid(True)

plt.figure()
plt.plot(x, uy.ravel())
plt.title("Uy along building line")
plt.xlabel("x [m]"); plt.ylabel("Uy [m]"); plt.grid(True)

plt.figure()
plt.plot(x, uz.ravel())
plt.title("Uz along building line")
plt.xlabel("x [m]"); plt.ylabel("Uz [m]"); plt.grid(True)

plt.show()

vertical_displacement_ground_building   = uz
horizontal_displacement_ground_building = ux

# Green-field analysis for Limiting Tensile Strain method:
print('--------- TENSILE STRAIN OUTPUT INFORMATION ---------')
dataSA = strain_analysis_greenfield(vertical_displacement_ground_building, horizontal_displacement_ground_building,
                                    length_beam_element, length_beam, building_height, neutral_line, poissons_ratio,
                                    input_type, def_mode)

print("Tensile strain calculation for Greenfield analysis")
highest_damage_greenfield, max_tensile_eps_gf = categorize_damage(dataSA['eps_t_max'])  # Highest damage category