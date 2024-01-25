"""
Objective: This code plots the U-data of magneticRotation 
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys

file_name = 'results.txt'
path = './'
# -------------------------- Comments ------------------------------------------
# input parameters specific file
# path to the data folder is ../data/data002_vd_20/files/ke_1024.txt for vd=20
# path to the data folder is ../data/data001_vd_80/files/ke_1024.txt for vd=80
# vd is an input parameter to run this file and that must be provided along with others.
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    t,ux,uy,uz = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([60,100]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot the figure
fig,ax = plt.subplots(2,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax[0].scatter(ux, uy, s=1.0)
ax[0].set_xlabel('ux')
ax[0].set_ylabel('uy')
ax[0].grid(True)

ax[1].plot(t,ux, linewidth=1.0, label='$u_{x}$')
ax[1].plot(t,uy, linewidth=1.0, label='$u_{y}$')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('ux, uy')
ax[1].legend(loc=0)
ax[1].grid(True)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.savefig(pjoin(path,'u.png'),dpi=dpi)
plt.show()
