import numpy as np 
import matplotlib.pyplot as plt 
from scipy.constants import value as constants 

e = constants('elementary charge')
me = constants('electron mass')
c = constants('speed of light in vacuum')

U = np.array([1E1, 1E2, 1E3, 1E4, 5E4, 1E5, 2E5, 1E6])
u = np.sqrt(2*e*U/me)
umax = np.max(u)

fig, ax = plt.subplots()
ax.semilogx(U, u/c, 'rd', linestyle='--', linewidth=0.5)
ax.set_xlabel('Potential (U)')
ax.set_ylabel('Final Velocity (u)')
ax.grid(True)
plt.show()

