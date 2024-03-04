import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys


file_name = 'ke.txt'
path = '../data/files/'

#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    t,keA,keB = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# -------- NORMALIZE THE INDIVIDUAL KINETIC ENERGIES W.R.T THE CODE------------
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
kB = constants('Boltzmann constant')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TA = keA/kB
TB = keB/kB

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([100,50]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 300                        #Print resolution
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
fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
hline = 750

ax.axhline(hline,linestyle='--',color='grey',linewidth=2.0)

ax.plot(t/1e-6,TA,'r',linestyle='-',linewidth=1.0,label='$T_{A}$')
ax.plot(t/1e-6,TB,'b',linestyle='-',linewidth=1.0,label='$T_{B}$')

ax.set_ylim(400, 1200)
ax.set_xlabel('t [$\mu$s]')
ax.set_ylabel('T [K]')
ax.legend(loc='upper right',framealpha=0.5)
ax.grid(True)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.savefig(pjoin('./temp.jpg'),dpi=dpi)

plt.show()
