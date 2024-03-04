"""
# This file is for the continuous time evolution of particle data.
# Only citing the path to the data file is enough to plot,
# No specific input corresponding to parameteric variation is required
# Run as: 
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin

path = '../data/'
path_fig = './'
# ------------------- Comments -------------------------------------------------
# Path to data file for vd=20 is
# ../data/particle_data/part_vd_20/data
# Path to data file for vd=80 is
# ../data/particle_data/part_vd_80/data
# ------------------------------------------------------------------------------
NUM_TS = 10000
write_interval = 100
rang = NUM_TS/write_interval

k = 0

for i in range(int(rang)):

    file_name_A = 'A%d.txt'%(int(k))
    file_name_B = 'B%d.txt'%(int(k))

    # Load A-Type particle data
    if os.path.exists(pjoin(path,file_name_A)):
        xA,vAx, vAy, vAz = np.loadtxt(pjoin(path,file_name_A),unpack=True)
    else:
        print('No data')
        exit()
    vA = np.sqrt(vAx**2 + vAy**2 + vAz**2)

    # Load B-Type particle data
    if os.path.exists(pjoin(path,file_name_B)):
        xB,vBx, vBy, vBz = np.loadtxt(pjoin(path,file_name_B),unpack=True)
    else:
        print('No data')
        exit()
    
    vB = np.sqrt(vBx**2 + vBy**2 + vBz**2)

    #---------------------------------------------------------
    # Plot the figure
    plt.cla()    
   
    plt.scatter(xA,vA,color='r',marker='.',s=1,alpha=0.5)
    plt.scatter(xB,vB,color='b',marker='.',s=1,alpha=0.5)    
    
    plt.xlabel('$x_{A,B}$')
    plt.ylabel('$\u03C5_{A,B}$')
    plt.xlim(0, 15E-6)
    plt.ylim(-1000, 2000)

    plt.title(i)
    print('k = ',k)
    k = k + write_interval
    #plt.savefig("./figure_part/part_%d.png"%(i))
    plt.pause(0.1)

plt.show()
