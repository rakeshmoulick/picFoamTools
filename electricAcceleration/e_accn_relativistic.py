import numpy as np 
import matplotlib as mpp
import matplotlib.pyplot as plt 
from scipy.constants import value as constants 
import os.path
from os.path import join as pjoin
from sympy import * 
# ------------------------------------------------------------------
e = constants('elementary charge')
me = constants('electron mass')
c = constants('speed of light in vacuum')
AMU = constants('atomic mass constant')
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
# ------------------------------------------------------------------
path = './'
# ------------------------------------------------------------------
def symbolic_coeffs():
    u, c, e, U, me = symbols('u, c, e, U, me')
    A = (2*e*U/me)**2       
    expr = (c**2)*(u**4) + A*(u**2) - A*(c**2)      
       
    p = Poly(expr, u)
    coff = p.all_coeffs()    
    return coff
# ------------------------------------------------------------------
coff = symbolic_coeffs()
# print(coff)
# print('Highest power of omega is: %d\n'%(len(coff)-1))
# ------------------------------------------------------------------
U = np.linspace(1E1, 1E6, 10000)
# U = np.array([1E1, 1E2, 1E3, 1E4, 5E4, 1E5, 2E5, 1E6])
# print(type(U))

def find_roots():
    # evaluate the coefficients which were symbolically found from the expression
    coeff1 = eval(str(coff[0]))
    coeff2 = eval(str(coff[1]))
    coeff3 = eval(str(coff[2]))
    coeff4 = eval(str(coff[3]))
    coeff5 = eval(str(coff[4]))
    
    # create an empty list of roots
    roots = []
    # form a vector of coefficients and then find root using np.roots(coeffs)
    for i in range(len(U)):
        coeffs = [coeff1, coeff2, coeff3[i], coeff4, coeff5[i]]
        root = np.roots(coeffs)
        # for each U there will be 6 roots, insert them serially in the empty array
        roots.append(root)
    
    # create a numpy array of the list of roots 
    roots = np.array(roots) 
    return roots

roots = find_roots()
# Find all the real roots
rr1 = np.real(roots[:,0])
rr2 = np.real(roots[:,1])
rr3 = np.real(roots[:,2])
rr4 = np.real(roots[:,3]) 
# Find all the imaginary roots
ri1 = np.imag(roots[:,0])
ri2 = np.imag(roots[:,1])
ri3 = np.imag(roots[:,2])
ri4 = np.imag(roots[:,3]) 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Data from Simulation (picFoam) -- Boris Particle Pusher
Usim_boris = np.array([1E1, 1E2, 1E3, 1E4, 5E4, 1E5, 2E5, 1E6])
usim_boris = np.array([1875541.827, 5930200.694, 18728216.08, 58456183.93, 123722016.5, 164354562.3, 208452051.3, 282128922.6])
# -----------------------------------------------------------------------------
# Data from Simulation (picFoam) -- Vay Particle Pusher
Usim_vay = np.array([1E1, 1E2, 1E3, 1E4, 5E4, 1E5, 2E5, 1E6])
usim_vay = np.array([1875541.827, 5930200.694, 18728216.08, 58456183.93, 123722016.5, 164354562.3, 208452051.3, 282128922.6])
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
figsize = np.array([100, 100/1.618]) # Figure size in mm
dpi = 300                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mpp.rc('text', usetex=False)
mpp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mpp.rc('axes', titlesize=10)
mpp.rc('axes', labelsize=10)
mpp.rc('xtick', labelsize=10)
mpp.rc('ytick', labelsize=10)
mpp.rc('legend', fontsize=10)
# ---------------------------------------------------------------------------
fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# we chose rr4 as this is the only real root. All other roots have imaginary parts
ax.semilogx(U, rr4/c, 'k', linestyle='--', linewidth = 0.5, label='Theory')
ax.semilogx(Usim_boris, usim_boris/c, 'rx', label='Simulation (Boris)')
ax.semilogx(Usim_vay, usim_vay/c, 'bd', label='Simulation (Vay)')

ax.set_xlabel('Potential (U) [V]')
ax.set_ylabel('Final Velocity (u/$c_{0}$) [-]')
ax.set_ylim([0, 1.2])
ax.set_xlim([1E1, 1E6])

ax.grid(True)
ax.legend(loc=0)
plt.savefig(pjoin(path,'electric_acceleration.png'),dpi=dpi)
plt.show()