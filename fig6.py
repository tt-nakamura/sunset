# reference:
#  A.T.Young, "Green Flashes and Mirages"
#   Optics and Photonics News, March 1999, p34

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos
from mirage import mirage

arcmin = np.pi/180/60

h,wavlen = 9.28,0.7
z = [0, 0.8, 1.8, 2.7, 4, 5.9, 8.9, 18.1, 32.5, 50, 1.1e4]
T = [288, 287.5, 287.3, 287.2, 287.1, 287, 286.9, 286.7, 286.5, 286.3, 216]

z_sun = [-15, -18, -21,
          -24, -27, -30,
          -33, -36, -39,
          -42, -45, -48,
          -50, -51, -52]
r_sun = 15.9*arcmin

A = atmos(z,T)
M = mirage(A, h, wavlen)

alpha = np.linspace(M.alpha_min, 30*arcmin, 500)
beta = alpha - M.refrac(alpha)
ymin = M.alpha_min/arcmin

plt.figure(figsize=(5,7.3))

for i,z in enumerate(z_sun):
    plt.subplot(5,3,i+1)
    plt.axis('image')
    plt.ylim([ymin,30])
    plt.xlim([-20,20])
    if i%3:  plt.tick_params(labelleft='off')
    if i<12: plt.tick_params(labelbottom='off')
    dz = beta - z*arcmin
    dy = np.where(np.abs(dz) > r_sun, np.nan,
                  (r_sun - dz)*(r_sun + dz))
    a = np.sqrt(dy)/arcmin

    plt.fill_betweenx(alpha/arcmin, -a, a,
                      color='red')
    plt.text(-19, 25, '%d'%z)

plt.tight_layout()
plt.savefig('fig6.png')
plt.show()
