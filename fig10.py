# reference:
#  A.T.Young, "Green Flashes and Mirages"
#   Optics and Photonics News, March 1999, p35

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos
from mirage import mirage

arcmin = np.pi/180/60

h,wavlen = 54,0.7
z = [0, 30, 50, 1e3, 1.1e4]
T = [287.3, 287.1, 287.8, 282, 215]

z_sun = [-16, -20, -24,
         -28, -32, -36,
         -40, -44, -48,
         -52, -56, -60,
         -62, -64, -65]
r_sun = 15.9*arcmin

A = atmos(z,T)
M = mirage(A, h, wavlen)

alpha = np.linspace(M.alpha_min, 30*arcmin, 500)
beta = alpha - M.refrac(alpha)
ymin = M.alpha_min/arcmin

plt.figure(figsize=(5,8.5))

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
plt.savefig('fig10.png')
plt.show()
