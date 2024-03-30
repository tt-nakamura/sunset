# references:
#  R.M.Green, "Spherical Astronomy" section 4.3
#  W.D.Bruton and G.W.Kattawar, "Unique temperature profile
#    for the atomosphere below an observer from sunset images"
#     Applied Optics 36 (1997) 6957

import numpy as np

class mirage:
    """ refraction of light in spherical earth """
    def __init__(self, atmos, height,
                 wavlen=0.6, rE=6378140, N=1024):
        """
        atmos: atoms object (see atmos.py)
        height: observer's height from ground in meter
        wavlen: wavelength in micro meter
        rE: earth radius in meter
        N: number of interpolation points for (nr, dlnn_dnr)
        assume nr increases with r
        """
        # R.M.Green, eq.(4.11)
        eps = 2.871e-4*(1 + 5.67e-3/wavlen**2)/1.29223
        self.atm = atmos
        self.eps = eps

        z = np.linspace(0, height, N)
        n, dlnn_dr = self.refrac_index(z)
        nr = n*(rE + z)
        dlnn_dnr = 1/(nr + n/dlnn_dr)
        a = -np.arccos(nr[0]/nr[-1])

        self.rE = rE
        self.h = height
        self.z = z
        self.n = n
        self.nr = nr
        self.dlnn_dnr = dlnn_dnr
        self.alpha_min = a # grazing ray with ground

    def refrac_index(self, z, dlnn_dz=True):
        """
        z: height from ground in meter
        dlnn_dz: bool, whether to ouput gradient of n
        if not dlnn_dz: return n = refractive index
        else: return n, dln(n)/dz
        z can be 1d-array (vectorized)
        """
        n1 = self.eps*self.atm.density(z) # n-1
        n = 1 + n1
        if not dlnn_dz: return n
        T = self.atm.temperature(z)
        dT = self.atm.lapse_rate(z) # dT/dz
        dlnn_dz = -n1/n/T*(self.atm.MgR + dT)
        return n, dlnn_dz

    def refrac_below(self, mu, N=256):
        """ refraction angle below the horizon
        mu: cos(alpha), alpha<0 (look down angle)
        N: number of trapezoids for integration
        return refraction angle in radian
        mu can be 1d-array (vectorized)
        """
        mu = np.asarray(mu)
        I = np.zeros_like(mu)
        if np.all(mu>=1): return I
        nr0 = self.nr[-1]
        nr1 = nr0*mu[mu<1]

        t = [np.linspace(0, np.sqrt(nr0-x), N) for x in nr1]
        t = np.stack(t, axis=1)
        nr = nr1 + t**2
        u = np.interp(nr, self.nr, self.dlnn_dnr)
        u = np.trapz(u/np.sqrt(nr + nr1), t, axis=0)
        I[mu<1] = -4*nr1*u
        return I

    def refrac_above(self, mu, N=64):
        """ refraction angle above the horizon
        mu: cos(alpha), alpha>0 (look up angle)
        N: number of trapezoids for integration
        return refraction angle in radian
        mu can be 1d-array (vectorized)
        """
        nr1 = self.nr[-1]*mu
        nr1_ = np.expand_dims(nr1, -1)
        def fun(z):
            n,dn = self.refrac_index(z)
            nr = n*(self.rE + z)
            return dn/np.sqrt((nr + nr1_)*(nr - nr1_))

        z_ = self.atm.z # atmospheric layers
        z_ = z_[z_>self.h]
        if len(z_)==0:
            raise RuntimeError(
                "oberver is too high or"
                "atmosphere is too low")

        # first layer above observer
        dt = np.sqrt(z_[0] - self.h)/N
        t = dt*(np.arange(N) + 0.5)
        z = self.h + t**2
        I = 2*np.sum(t*fun(z), axis=-1)*dt

        # intermediate layers
        for za,zb in zip(z_[:-1], z_[1:]):
            z = np.linspace(za,zb,N)
            I += np.trapz(fun(z), z)

        # top layer (exponentially decreasing density)
        t = np.arange(N) + 0.5
        z = z_[-1]*(1 - np.log(t/N))
        I += np.sum(fun(z)/t, axis=-1)*z_[-1]
        I *= -nr1

        return I

    def refrac(self, alpha, N=64, M=256):
        """ refraction angle
        alpha: elevation angle / radian
        N: number of trapezoids (passed to refrac_above)
        M: number of trapezoids (passed to refrac_below)
        return refraction angle in radian
        alpha can be 1d-array (vectorized)
        assume alpha > alpha_min (see __init__())
        """
        a = np.asarray(alpha)
        if np.any(a < self.alpha_min):
            raise RuntimeError("negative alpha too large")
        
        mu = np.cos(a)
        I = self.refrac_above(mu, N)
        I = np.asarray(I)

        I[a<0] += self.refrac_below(mu[a<0], M)
        return I

    def AbelInv(self, data, z, N=256):
        """ infer refractive index from refaction data
            using Abel's integral inversion theorem
        data: (alpha, F) tuple, where
          alpha: 1d-array, elevation angle in radian
          F: 1d-array, output of rafrac_below
          assume alpha is negative and increasing
          assume len(alpha) == len(F)
        z: height at which to infer temperature
        N: number of trapezoids for integration
        return refractive index at z
        z can be 1d-array (vectorized)
        if z is 1d-array, assume z increases
        assume z < observer's height
        """
        if np.any(z > self.h):
            raise RuntimeError("z>h in AbelInv")

        alpha,F = data
        mu = np.cos(alpha)
        if np.any(np.diff(mu) < 0):
            raise RutimeError("mu must increase")
            
        if np.isscalar(z):
            return self.AbelInv_(mu, F, z, None, N)

        n,n0 = [],None
        for z in z[::-1]:
            n0 = self.AbelInv_(mu, F, z, n0, N)
            n.append(n0)
        return np.asarray(n[::-1])

    def AbelInv_(self, mu, F, z, n_init, N):
        """ private function used in AbelInv
        mu: cos(alpha)
        F: output of refrac_below
        z: oberver's height, scalar
        n_init: intial guess for refractive index
        N: number of trapezoids for integration
        return refractive index at z
        ref. Bruton and Kattawar, eq.(13)
        """
        n0 = self.n[-1]
        if z == self.h: return n0
        nr0 = self.nr[-1]
        n = n_init if n_init else n0 # initial guess
        r = z + self.rE

        while True:
            nr = n*r
            t = np.linspace(0, np.sqrt(nr0-nr), N)
            nr1 = nr + t**2
            u = np.interp(nr1/nr0, mu, F)
            I = np.trapz(u/np.sqrt(nr1 + nr), t)
            s = n0*np.exp(2/np.pi*I)
            if np.isclose(s-1,n-1): return s
            n = (n+s)/2 # moderate update

    def RhoFromN(self, n):
        """ density from refractive index """
        return (n-1)/self.eps # kg/m^3
