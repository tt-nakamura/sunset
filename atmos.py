# reference:
#  "Manned Spacecraft Design Principles" chapter 2

import numpy as np
from scipy.constants import g,R
from scipy.integrate import cumtrapz

class atmos:
    """ earth atmosphere in gravitational equilibrium """
    def __init__(self, z, T=None, dTdz=None,
                 T0=273.15, p0=1.01325e5, M=28.964e-3):
        """
        z: height from ground in meter (1d-array)
        T: temperature in K (1d-array)
        dTdz: temperature gradient in K/m (1d-array)
        T0: temperature at z[0] in K
        p0: pressure at z[0] in Pascal
        M: mean molecular weight of air in kg/mol
        if T is not None, dTdz and T0 are ignored
        otherwise, dTdz and T0 must be supplied
        assume len(T)==len(z) or len(dTdz)==len(z)-1
        assume z increases
        isothermal at z higher than z[-1]
        """
        z = np.asarray(z)
        dz = np.diff(z)
        if T is not None:
            if len(T) != len(z):
                raise RuntimeError("bad len(T)")
            T = np.asarray(T)
            dT = np.r_[np.diff(T)/dz, 0]
        elif dTdz is not None:
            if len(dTdz) != len(dz):
                raise RuntimeError("bad len(dTdz)")
            T = np.cumsum(np.r_[T0, dz*dTdz])
            dT = np.r_[dTdz, 0]
        else:
            raise RuntimeError("T and dTdz are None")

        MgR = M*g/R
        MgRT = MgR/T
        MgRdT = MgR/np.where(dT, dT, 1)
        p = (T[:-1]/T[1:])**MgRdT[:-1]
        p_ = np.exp((z[:-1] - z[1:])*MgRT[:-1])
        p = np.where(dT[:-1], p, p_)
        p = np.cumprod(np.r_[p0, p])
        rho = p*M/R/T

        self.dT = dT
        self.z = z
        self.T = T
        self.p = p
        self.rho = rho
        self.MgR = MgR
        self.MgRT = MgRT
        self.MgRdT = MgRdT

    def temperature(self, z):
        """ temperature at z
        z: height in meter (scalar or 1d-array)
        return temperature in K
        """ # linear interpolation
        return np.interp(z, self.z, self.T)

    def lapse_rate(self, z):
        """ temperature gradient dT/dz at z
        z: height in meter (scalar or 1d-array)
        return dT/dz in K/m
        """
        i = np.searchsorted(self.z, z, 'r') - 1
        return self.dT[i]

    def pressure(self, z):
        """ pressure at z
        z: height in meter (scalar or 1d-array)
        return pressure in Pascal
        """
        i = np.searchsorted(self.z, z, 'r') - 1
        t = self.temperature(z)/self.T[i]
        t = t**(-self.MgRdT[i])
        t_ = np.exp(-(z - self.z[i])*self.MgRT[i])
        t = np.where(self.dT[i], t, t_)
        return self.p[i]*t

    def density(self, z):
        """ density at z
        z: height in meter (scalar or 1d-array)
        return density in kg/m^3
        """
        i = np.searchsorted(self.z, z, 'r') - 1
        t = self.temperature(z)/self.T[i]
        t = t**(-self.MgRdT[i] - 1)
        t_ = np.exp(-(z - self.z[i])*self.MgRT[i])
        t = np.where(self.dT[i], t, t_)
        return self.rho[i]*t
   
    def TFromRho(self, rho, z):
        """ compute T(z) from rho(z)
        z: height in meter (1d-array)
        rho: denisty in kg/m^3 (1d-array)
        assume len(z)==len(rho)
        """ # integrate dp/dz = -g*rho numerically
        p = self.p[0]/g - cumtrapz(rho, z, initial=0)
        return p/rho*self.MgR
