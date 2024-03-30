# reference:
#  A.T.Young, "Green Flashes and Mirages"
#   Optics and Photonics News, March 1999, p33

import numpy as np
import matplotlib.pyplot as plt
from atmos import atmos

z = [0, 0.8, 1.8, 2.7, 4, 5.9, 8.9, 18.1, 32.5, 50, 1.1e4]
T = [288, 287.5, 287.3, 287.2, 287.1, 287, 286.9, 286.7, 286.5, 286.3, 215]

A = atmos(z,T)

z = np.linspace(0,100,100)

plt.figure(figsize=(5,5))

plt.subplot(311)
plt.plot(z, A.temperature(z))
plt.ylabel('temperature / K')
plt.tick_params(labelbottom='off')

plt.subplot(312)
plt.plot(z, A.pressure(z)/1e3)
plt.ylabel('pressure / kPa')
plt.tick_params(labelbottom='off')

plt.subplot(313)
plt.plot(z, A.density(z))
plt.ylabel('density / kg/m$^3$')
plt.xlabel('$z$ = height / m')

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
