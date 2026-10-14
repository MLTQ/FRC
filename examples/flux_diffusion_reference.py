"""Independent continuum solution for a prescribed sin-squared enclosed flux."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.special import j0, j1, jn_zeros


class FluxDiffusionReference:
    def __init__(self, radial_points, *, radius, density, ion_mass, resistivity,
                 amplitude, pulse_duration, modes=512):
        self.r = np.asarray(radial_points, float)
        self.radius = radius; self.amplitude = amplitude; self.period = pulse_duration
        self.omega = 2*np.pi/pulse_duration
        reduced = ME*ion_mass/(ME+ion_mass)
        self.de2 = reduced/(MU0*density*E**2)
        self.diffusivity = resistivity/MU0
        roots = jn_zeros(1, modes); self.q = roots/radius
        self.jwall = j0(roots)
        self.rates = self.diffusivity*self.q**2/(1+self.de2*self.q**2)
        self.basis_b = j0(self.q[:, None]*self.r.ravel()[None, :])
        self.basis_j = self.q[:, None]*j1(self.q[:, None]*self.r.ravel()[None, :])/MU0

    def convolution(self, rates, when):
        """Integral exp[-lambda(t-s)] Bbar_dot(s) ds, evaluated analytically."""
        rates = np.asarray(rates)
        t = min(max(float(when), 0.), self.period)
        if t == self.period:
            sine = 0.; one_minus_cosine = 0.
        else:
            sine = np.sin(self.omega*t)
            one_minus_cosine = 2*np.sin(.5*self.omega*t)**2
        value = self.amplitude*np.pi/self.period*(rates*sine
            +self.omega*(np.expm1(-rates*t)+one_minus_cosine))/(rates**2+self.omega**2)
        return value*np.exp(-rates*max(float(when)-self.period, 0.))

    def evaluate(self, when):
        t = float(when)
        mean = self.amplitude*np.sin(np.pi*t/self.period)**2 if 0 < t < self.period else 0.
        # Exact polynomial lifting retains the nonzero physical wall current;
        # the residual series converges rapidly instead of approximating its jump.
        a = float(self.convolution(self.diffusivity/self.de2, t)/(4*self.de2))
        modal = self.convolution(self.rates, t)/(self.jwall*(1+self.de2*self.q**2))
        residual = modal-4*a/(self.q**2*self.jwall)
        bz = mean+a*(self.r**2-.5*self.radius**2)+(residual@self.basis_b).reshape(self.r.shape)
        jtheta = -2*a*self.r/MU0+(residual@self.basis_j).reshape(self.r.shape)
        return bz, jtheta
