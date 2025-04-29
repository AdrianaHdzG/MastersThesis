from sys import platform as pltm
import os, platform, pkg_resources
from ctypes import CDLL, c_int, c_double, c_char_p, POINTER
import numpy as np
import _ctypes
import ASREpy.ground_deformation_models as gdm
import matplotlib.pyplot as plt

import importlib.util
import sys
import os

module_name = "ASRE_Timoshenko_model"
module_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "ASRE_Timoshenko_model.py")

spec = importlib.util.spec_from_file_location(module_name, module_path)
module = importlib.util.module_from_spec(spec)
sys.modules[module_name] = module
spec.loader.exec_module(module)

# class ASRE_Timoshenko_model:
#     """
#     A class to represent an ASRE Timoshenko beam model.
    
#     Attributes
#     ----------
#     nnode : int
#         Number of node.
#     meshX : np.array(dtype=np.float64)
#         The x coordinates of the beam nodes. Unit: m
#     meshY : np.array(dtype=np.float64)
#         The y coordinates of the beam nodes. Unit: m
#     meshZ : np.array(dtype=np.float64)
#         The z coordinates of the beam nodes. Unit: m
#     dfoot : float
#         The depth of the beam. Unit: m
#     bfoot : float
#         The width of the beam. Unit: m
#     Eb : float
#         The Young's modulus of the beam. Unit: N/m^2
#     EoverG : float
#         The ratio between Young's modulus and the shear modulus. Unit: unitless
#     ni_foot : float
#         The Poisson's ratio of beam. Unit: unitless
#     q_foot : float
#         The uniform weight applied on the beam. Unit: N/m (Weight along unit length in the longitudinal direction)
#     EsNominal : float
#         The nominal elastic modulus of soil. Unit: N/m^2
#     nis : float
#         The Poisson's ratio of soil. Unit: unitless 
#     mu_int : float
#         The friction coefficient between soil and beam. Unit: unitless
#     """

#     def __init__(self, nnode, meshX, meshY, meshZ, dfoot, bfoot):
#         """
#         Constructs the beam with dimension properties.

#         Parameters
#         ----------
#         nnode : int
#             Number of node.
#         meshX : np.array(dtype=np.float64)
#             The x coordinates of the beam nodes. Unit: m
#         meshY : np.array(dtype=np.float64)
#             The y coordinates of the beam nodes. Unit: m
#         meshZ : np.array(dtype=np.float64)
#             The z coordinates of the beam nodes. Unit: m
#         dfoot : float
#             The depth of the beam. Unit: m
#         bfoot : float
#             The width of the beam. Unit: m
#         """
#         self.nnode = nnode
#         self.meshX = meshX
#         self.meshY = meshY
#         self.meshZ = meshZ
#         self.dfoot = dfoot
#         self.bfoot = bfoot
#         self.asre_dll = self._import_dll()

#     def _import_dll(self):
#         """
#         import the CDLL.
        
#         Returns
#         -------
#         CDLL
#             The CDLL.
#         """
#         if pltm == "linux" or pltm == "linux2":
#             c_lib = CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                                     "lib/linux/ASRElib.so"))
#         elif pltm == "darwin":
#             if platform.processor() == 'arm':
#                 c_lib = CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                                         "lib/macOS_m1/ASRElib.so"))
#             else:
#                 c_lib = CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                                         "lib/macOS/ASRElib.so"))
#         elif pltm == "win32":
#             c_lib = CDLL(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
#                                     "bin", "win32", "ASRElib.dll"))
#             # c_lib = CDLL(pkg_resources.resource_filename('ASREpy', 'ASREcpp//bin//win32//ASRElib.dll'))
#             # c_lib.printName()
#         c_lib.run.argtypes = [c_int, #nnode 
#                               np.ctypeslib.ndpointer(dtype=np.float64), #meshX
#                               np.ctypeslib.ndpointer(dtype=np.float64), #meshY
#                               np.ctypeslib.ndpointer(dtype=np.float64), #meshZ
#                               np.ctypeslib.ndpointer(dtype=np.float64), #dispV
#                               np.ctypeslib.ndpointer(dtype=np.float64), #dispL
#                               np.ctypeslib.ndpointer(dtype=np.float64), #dispT
#                               c_double, #Eb
#                               c_double, #EoverG
#                               c_double, #EsNominal
#                               c_double, #nis
#                               c_double, #dfoot
#                               c_double, #bfoot
#                               c_double, #ni_foot
#                               c_double, #mu_int
#                               c_double, #qz_foot
#                               c_char_p, #output
#                                 #   POINTER(c_double)
#                                 #   np.ctypeslib.ndpointer(dtype=np.float64)
#                             ]
#         return c_lib

#     def set_beam_properties(self, Eb, EoverG, q_foot):
#         """
#         Set the beam material properties.

#         Parameters
#         ----------
#         Eb : float
#             The Young's modulus of the beam. Unit: N/m^2
#         EoverG : float
#             The ratio between Young's modulus and the shear modulus. Unit: unitless
#         q_foot : float
#             The uniform weight applied on the beam. Unit: N/m (Weight along unit
#             length in the longitudinal direction)
#         """
#         self.Eb = Eb
#         self.EoverG = EoverG
#         self.q_foot = q_foot
#         self.ni_foot = self.EoverG/2 - 1 #The Poisson's ratio of beam. Unit: unitless

#     def set_soil_properties(self, EsNominal, nis, mu_int):
#         """
#         Set the beam material properties.

#         Parameters
#         ----------
#         EsNominal : float
#             The nominal elastic modulus of soil. Unit: N/m^2
#         nis : float
#             The Poisson's ratio of soil. Unit: unitless 
#         mu_int : float
#             The friction coefficient between soil and beam. Unit: unitless
#         """
#         self.EsNominal = EsNominal
#         self.nis = nis
#         self.mu_int = mu_int
        
#     def run_model(self, dispL, dispT, dispV, output = 'disp'):
#         """
#         Run the SSI model under greenfield displacements

#         Parameters
#         ----------
#         dispL : float
#             The nominal elastic modulus of soil. Unit: N/m^2
#         dispT : float
#             The Poisson's ratio of soil. Unit: unitless 
#         dispV : float
#             The friction coefficient between soil and beam. Unit: unitless
#         output : str
#             If 'disp' then save beam displacement to self.beam_disp
#             self.beam_dispL is the beam disp in longitudinal direction
#             self.beam_dispT is the beam disp in transverse direction
#             self.beam_dispV is the beam disp in vertical direction
#             If 'strain' then save beam principal strain to self.strain
#             Each element in self.strain is one principal strain of the beam
#             If 'strain+disp' then save both self.disp and self.strain
        
#         Returns
#         -------
#         bool
#             Return Ture if run success and False if unsuccess 
#         """
#         self.ouput = output
#         if self.ouput == 'disp':
#             self.asre_dll.run.restype = POINTER(c_double * (self.nnode * 6))
#         elif self.ouput == 'strain':
#             self.asre_dll.run.restype = POINTER(c_double * 3)
#         elif self.ouput == 'strain+disp':
#             self.asre_dll.run.restype = POINTER(c_double * (3 + self.nnode * 6))
#         else:
#             raise ValueError(f'output value {output} is not permitted in ASRE_Timoshenko_model')
        
#         self.ouput = self.ouput.encode('utf-8')
#         try:
#             result = self.asre_dll.run(self.nnode, self.meshX, self.meshY, self.meshZ,
#                                        dispV, dispL, dispT, self.Eb, self.EoverG,
#                                        self.EsNominal, self.nis, self.dfoot,
#                                        self.bfoot, self.ni_foot, self.mu_int,
#                                        self.q_foot, self.ouput)
#             result_list = [i for i in result.contents]
#             # print(type(result_list))
#             # print(result_list)
#             result_array = np.array(result_list)
#             # print(type(result_array))
#             # print(result_array)
#             if self.ouput.decode('utf-8') == 'disp':
#                 self.beam_DispL = result_array[0::6]
#                 self.beam_DispT = result_array[1::6]
#                 self.beam_DispV = result_array[2::6]
#                 self.beam_RotaL = result_array[3::6]
#                 self.beam_RotaT = result_array[4::6]
#                 self.beam_RotaV = result_array[5::6]
#             elif self.ouput.decode('utf-8') == 'strain':
#                 self.beam_p_strains = result_array
#             elif self.ouput.decode('utf-8') == 'strain+disp':
#                 self.beam_p_strains = result_array[0:3]
#                 result_array = result_array[3:]
#                 self.beam_DispL = result_array[0::6]
#                 self.beam_DispT = result_array[1::6]
#                 self.beam_DispV = result_array[2::6]
#                 self.beam_RotaL = result_array[3::6]
#                 self.beam_RotaT = result_array[4::6]
#                 self.beam_RotaV = result_array[5::6]
#             return True
#         except:
#             self.release_cdll_handle()
#             raise RuntimeError(f'ASRE_Timoshenko_model failed to run the ASRE cpp library')
    
#     def release_cdll_handle(self):
#         if self.asre_dll is None:
#             pass
#         else:
#             lib_handle = self.asre_dll._handle
#             del self.asre_dll
#             _ctypes.FreeLibrary(lib_handle)
#             self.asre_dll = None

 
if __name__=="__main__":
    # STR-1, vl = 0.5
    beam_id = "STR-1"
    beamX = np.linspace(-15, 15, 61)
    z0 = 11.25
    vl = 0.5/100
    k = 0.55
    D = 6.16
    dispV, dispH = gdm.ground_disp_Mair_1993(beamX, z0, k, D, vl)
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True,
                                    figsize=(12, 6))
    axes[0].plot(beamX, dispH, label = 'Gaussian')
    axes[1].plot(beamX, dispV, label = 'Gaussian')
    dfoot = 0.12
    bfoot = 10
    beamY = np.zeros_like(beamX)
    beamZ = np.zeros_like(beamX)
    Eb = 70e9
    EoverG = 2.6 # Assume very large G to model the Euler–Bernoulli beam used in Franza and DeJong
    qfoot = 3.2*10*1000
    Es = 25e6
    nis = 0.25
    mu_int = np.tan(30*np.pi/180)
    model = module.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                        beamZ, dfoot, bfoot)
    model.set_beam_properties(Eb, EoverG, qfoot)
    model.set_soil_properties(Es, nis, mu_int)

    model.run_model(dispH, np.zeros_like(dispH), dispV, "disp")
    print(model.beam_DispV)
