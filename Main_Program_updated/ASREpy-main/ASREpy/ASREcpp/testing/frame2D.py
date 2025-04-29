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
    EoverG = 0.001 # Assume very large G to model the Euler–Bernoulli beam used in Franza and DeJong
    qfoot = 3.2*10*1000
    Es = 25e6
    nis = 0.25
    mu_int = np.tan(30*np.pi/180)
    model = module.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                        beamZ, dfoot, bfoot)
    model.set_beam_properties(Eb, EoverG, qfoot)
    model.set_soil_properties(Es, nis, mu_int)

    cdll = model.asre_dll

    cdll.calFlexVaziri.argtypes = [ np.ctypeslib.ndpointer(dtype=np.float64),
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    c_double,
                                    c_int,
                                    np.ctypeslib.ndpointer(dtype=np.float64),
                                    c_double,
                                    np.ctypeslib.ndpointer(dtype=np.float64)]
    # cdll.calFlexVaziri.restype = POINTER(c_double * model.nnode * model.nnode)

    cdll.calFlexVaziri.restype = POINTER(c_double * (2*3) * (2*3))
    FLEX = cdll.calFlexVaziri(model.meshZ[0:2], model.meshX[0:2], model.meshY[0:2], Es, 2,
                       np.ones_like(model.meshZ[0:2])*0.5, nis, np.ones_like(model.meshZ[0:2])*dfoot)
    
    FLEX_list = [i for i in FLEX.contents]
    print(np.ctypeslib.as_array(FLEX.contents, (6,6)))
    print(np.array(FLEX_list[0]))