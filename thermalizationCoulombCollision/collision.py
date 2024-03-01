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
# y(0) = e, y(1) = i
def F(y, t):
    F = np.zeros(2)
    F[0] = nu*(y[1] - y[0])
    F[1] = nu*(y[0] - y[1])
    return F
# -----------------------------------------
# def coll():
#     y0 = [Te, Ti]
#     t = np.linspace(0, 2E-6, 100)
#     y = odeint(F, y0, t)
#     return t, y
# -----------------------------------------
Te = 20 # electrons
Ti = 10 # Ar-ions
T_bar = (Te + Ti)/2.0
# -----------------------------------------
me = me
mi = 39.948*amu
amug = amu*1E3 # amu in grams
# -----------------------------------------
ne = 3.034E22
ni = 3.034E22
# -----------------------------------------
Ze = 1
Zi = 1
Z = 1
ln_lambda = 15
# -----------------------------------------
# A1 = (Ze**2)*(Zi**2)*(ni*1e6)*np.sqrt(me*1E3*mi*1E3)*ln_lambda
# A2 = (me*1E3*Te + mi*1E3*Ti)**(3/2)
# A1 = (Ze**2)*(Zi**2)*(ni*1E6)*np.sqrt(me/amug*mi/amu)*ln_lambda
# A2 = (me/amug*Te + mi/amu*Ti)**(3/2)
# nu = 1.8E-19*A1/A2

# nu = nu*1.5E3
# nu = 2.2e9
# --------------------------------------------------
vte = np.sqrt(Te*e/me)
A1 = pi**(3/2)*ne*Z*(e**4)*ln_lambda
A2 = np.sqrt(2)*((4*pi*eps0)**2)*(me**2)*(vte**3)
nu = A1/A2 
# --------------------------------------------------
print(nu)

def main():    
    # print(nu)
    # -----------------------------------------
    y0 = [Te, Ti]
    t = np.linspace(0, 2E-9, 100)
    y = odeint(F, y0, t)
    plt.plot(t/1e-9,y[:,0], label='$e$')
    plt.plot(t/1e-9,y[:,1], label='$Ar+$')
    plt.axhline(15, color='grey', linestyle='--')
    plt.xlabel('Time[ns]', fontsize = 15)
    plt.ylabel('Temperature[eV]', fontsize = 15)
    plt.xticks(fontsize = 15) 
    plt.yticks(fontsize = 15) 
    plt.legend(loc=0, fontsize=15)
    plt.grid(True)
    plt.savefig('thermalization_theory', dpi=300)
    plt.show()

if __name__== "__main__":
    main()
