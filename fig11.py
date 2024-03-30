# reference:
#  W.D.Bruton and G.W.Kattawar, "Unique temperature profile
#    for the atomosphere below an observer from sunset images"
#     Applied Optics 36 (1997) 6957

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos
from mirage import mirage

h = 54
z = [0, 30, 50, 1e3, 1.1e4]
T = [287.3, 287.1, 287.8, 282, 215]

A = atmos(z,T)
M = mirage(A, h, 0.7)

alpha = np.linspace(M.alpha_min, 0, 100)
F = M.refrac(alpha) - M.refrac(-alpha)

z = np.linspace(0, h, 100)
n = M.AbelInv((alpha,F), z)
rho = M.RhoFromN(n)
T = A.TFromRho(rho,z)

plt.figure(figsize=(5, 3.75))
plt.plot(z, A.temperature(z), label='fig5')
plt.plot(z[::3], T[::3],'.', label="Abel's theorem")
plt.xlabel('$z$ = height / m')
plt.ylabel('$T$ = temperature / K')
plt.legend()

plt.tight_layout()
plt.savefig('fig11.eps')
plt.show()
