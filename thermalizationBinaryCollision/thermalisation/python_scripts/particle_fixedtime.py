"""
# This file plots the particle data at fixed times. This is the file to generate plots for the paper.
# Only citing the path to the data file is enough to plot,
# No specific input corresponding to parameteric variation is required
# Run as: 
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
import seaborn as sns
sns.set(style='whitegrid')
import os.path
from os.path import join as pjoin
import sys

#path = sys.argv[1]
path = './data'
path_fig = './'
# ------------- Comments -------------------------------------------------------
# Path to data file for vd=20 is
# ../data/particle_data/part_vd_20/data
# Path to data file for vd=80 is
# ../data/particle_data/part_vd_80/data
# Figures to be saved in ../data/particle_data/part_vd_20/figs for vd=20 &
# ../data/particle_data/part_vd_80/figs for vd=80
#-------------------------------------------------------------------------------
# Define the time at which data is required

#write_interval = 200
# File index value is k(say), change this index to get the plot at required time.
# Calculate the wpet for a k-value by using DT_coeff beforehand.
k = 2500

file_name_A = 'A%d.txt'%(int(k))
file_name_B = 'B%d.txt'%(int(k))

# Load Cold electron data
if os.path.exists(pjoin(path,file_name_A)):
    xA,vxA, vyA, vzA = np.loadtxt(pjoin(path,file_name_A),unpack=True)
else:
    print('No data')
    exit()

# Load Hot electron data
if os.path.exists(pjoin(path,file_name_B)):
    xB,vxB, vyB, vzB = np.loadtxt(pjoin(path,file_name_B),unpack=True)
else:
    print('No data')
    exit()

# Plot the figure
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([150,200]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                          #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24   #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# In scatter plot s specifies the area
fig,ax = plt.subplots(2,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# --- Play with the parameter 's' in the plot. For vd=20 at wpet=100, 200 make s = 0.00001 for clarity else use 0.01/0.001

#plt.figure(1)
ax[0].scatter(xA,vzA,color='g',marker='.',edgecolor='grey',s=1.0)
ax[0].set_xlabel('$x_{A}$')
ax[0].set_ylabel('$\u03C5_{z}A$')

#plt.figure(2)
ax[1].scatter(xB,vxB,color='blue',marker='.',edgecolor='blue',s=1.0)
ax[1].set_xlabel('$x_{B}$')
ax[1].set_ylabel('$\u03C5_{z}B$')


# Save Figure
plt.savefig(pjoin(path_fig,'part.png'),format='png',dpi=300)

plt.show()
