import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
# import ASRE models
cur_dir = os.getcwd()
ASREpy_dir = os.path.join(os.path.dirname(os.path.dirname(cur_dir)))
sys.path.append(ASREpy_dir)

# For Sebastian
sys.path.append(r"C:\Users\sebas\OneDrive - Aarhus universitet\project_2025_sebastian_peter\03_Code\Main_Program\ASREpy-main")

# For Peter
#sys.path.append(r"C:\Users\Peter\OneDrive - Aarhus universitet\project_2025_sebastian_peter\03_Code\Main_Program\ASREpy-main")

import ASREpy
import ASREpy.ground_deformation_models as gdm
import pandas as pd
from scipy import interpolate
import scipy
np.set_printoptions(threshold=np.inf)