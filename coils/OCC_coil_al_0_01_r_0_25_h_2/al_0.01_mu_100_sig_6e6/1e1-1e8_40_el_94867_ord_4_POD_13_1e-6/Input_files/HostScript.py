import os
import numpy as np
#from contextlib import chdir

#with chdir(Your_Path):
# do stuff here
from main import main
import time

if __name__ == '__main__':
    comparison_eig_p = np.zeros((40,4), dtype=complex)
    comparison_ndofs_p = np.zeros(4)
    CPUs=[3,3,2,2,1]
    for p in [0,1,2,3,4]:
        print("order running = ",p)
        Return_Dict = main(geometry='OCC_coil_al_0_01_r_0_25_h_2.py',cpus=CPUs[p],use_OCC=True, use_POD=True, order=p)
