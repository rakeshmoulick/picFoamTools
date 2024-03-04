import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import odeint 
from scipy.constants import value as constants
import scipy
# -----------------------------------------
eps0 = constants('electric constant')
kB   = constants('Boltzmann constant')
me   = constants('electron mass')
amu  = constants('atomic mass constant')
e    = constants('elementary charge')
pi   = scipy.constants.pi
# -----------------------------------------
def F(y, t):
    F = np.zeros(2)
    F[0] = nu_ab*(y[1] - y[0])
    F[1] = nu_ba*(y[0] - y[1])
    return F

def coll():
    y0 = [Ta, Tb]
    t = np.linspace(0, 2E-6, 1000)
    y = odeint(F, y0, t)
    return t, y
# -----------------------------------------
Ta = 1000 # Ne
Tb = 500  # Ar
T_bar = (Ta + Tb)/2.0
# -----------------------------------------
ma = 20.18*amu
mb = 39.948*amu
# -----------------------------------------
na = 1E22
nb = 1E22
# -----------------------------------------
rad_a = 154*1E-12; rad_b = 188*1E-12; 
dia_a = 2*rad_a; dia_b = 2*rad_b
sigmaT = (pi/4)*(dia_a + dia_b)**2
print("sigma: ", sigmaT)
# -----------------------------------------
mu = ma*mb/(ma+mb)
gBar = np.sqrt((8*kB*T_bar)/(pi*mu))
# -----------------------------------------
nu_ab = nb*gBar*sigmaT
nu_ba = na*gBar*sigmaT
# -----------------------------------------

def main():
    t, y = coll()
    plt.plot(t/1E-6,y[:,0], label='$T_{\u03B1}$')
    plt.plot(t/1E-6,y[:,1], label='$T_{\u03B2}$')
    hline = T_bar
    plt.axhline(hline, color='grey', linestyle='--', label='$T_{avg}$')

    plt.xlabel('Time [$\mu$s]')
    plt.ylabel('Temperature [K]')
    plt.legend(loc=0)
    plt.grid(True)
    plt.show()

if __name__== "__main__":
    main()
