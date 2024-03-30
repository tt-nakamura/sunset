# reference:
#  A.T.Young, "Green Flashes and Mirages"
#   Optics and Photonics News, March 1999, p33

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos
from mirage import mirage

arcmin = np.pi/180/60

h = 9.28
z = [0, 0.8, 1.8, 2.7, 4, 5.9, 8.9, 18.1, 32.5, 50, 1.1e4]
T = [288, 287.5, 287.3, 287.2, 287.1, 287, 286.9, 286.7, 286.5, 286.3, 215]

A = atmos(z,T)
R = mirage(A, h, 0.7)
G = mirage(A, h, 0.5)

x = np.linspace(R.alpha_min, 40*arcmin, 200)
yr = R.refrac(x)
yg = G.refrac(x)

plt.figure(figsize=(5,3.75))
plt.plot(x/arcmin, yr/arcmin, 'red',
         label=r'$\lambda = 0.7\,\mu$m (red)')
plt.plot(x/arcmin, yg/arcmin, 'lime',
         label=r'$\lambda = 0.5\,\mu$m (green)')
plt.xlabel(r'$\alpha$ = apparent elevation / arcmin')
plt.ylabel(r'$F(\alpha)$ = angle of refraction / arcmin')

plt.legend()
plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
