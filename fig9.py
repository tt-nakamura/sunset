# reference:
#  A.T.Young, "Green Flashes and Mirages"
#   Optics and Photonics News, March 1999, p35

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos
from mirage import mirage

arcmin = np.pi/180/60

h = 54
z = [0, 30, 50, 1e3, 1.1e4]
T = [287.3, 287.1, 287.8, 282, 215]

A = atmos(z,T)
R = mirage(A, h, 0.7)
G = mirage(A, h, 0.5)

x = np.linspace(G.alpha_min, 40*arcmin, 200)
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
plt.savefig('fig9.eps')
plt.show()
