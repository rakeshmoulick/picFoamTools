import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys

file_name_1 = "vdf0.txt"
path = './'

if os.path.exists(pjoin(path,file_name_1)):
    v0, count0 = np.loadtxt(pjoin(path,file_name_1), unpack=True)    
else:
    print('No data')
    exit()


file_name_2 = "vdf1.txt"
path = './'

if os.path.exists(pjoin(path,file_name_2)):
    v, count = np.loadtxt(pjoin(path,file_name_2), unpack=True)    
else:
    print('No data')
    exit()


plt.plot(v0, count0)
plt.plot(v, count)
plt.grid(True)

plt.show()