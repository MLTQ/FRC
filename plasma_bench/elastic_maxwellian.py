"""Conservative finite-drift elastic moments from a momentum-transfer table."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import elementary_charge as EV
from scipy.stats import ncx2


class TabulatedElasticPair:
    """Two isotropic Maxwellians; only the supplied relative-speed domain enters.

    Table energy is the incident projectile energy in the stationary target
    frame, E=mass_a*g**2/2. ELASTIC must mean momentum-transfer cross section.
    No extrapolated cross section or assumed missing-tail contribution is added.
    """
    def __init__(self, section, mass_a, mass_b, *, order=16):
        if section.kind != 'ELASTIC':
            raise ValueError('An explicit elastic momentum-transfer table is required')
        if not np.isfinite([mass_a, mass_b]).all() or min(mass_a, mass_b) <= 0:
            raise ValueError('Positive finite physical masses required')
        if isinstance(order, bool) or not isinstance(order, int) or order < 2:
            raise ValueError('Integer quadrature order >=2 required')
        e, sigma = np.asarray(section.energy_ev), np.asarray(section.sigma_m2)
        if (e.ndim != 1 or e.shape != sigma.shape or len(e) < 2 or
                not np.isfinite(e).all() or not np.isfinite(sigma).all() or
                np.any(e < 0) or np.any(sigma < 0) or np.any(np.diff(e) <= 0)):
            raise ValueError('Finite nonnegative table with increasing energies required')
        self.section = section
        self.mass_a, self.mass_b = float(mass_a), float(mass_b)
        self.speed_edges = np.sqrt(2*EV*e/self.mass_a)
        self.nodes, self.weights = leggauss(order)
        self.order = order

    def coefficients(self, temperature_a_ev, temperature_b_ev, drift_speed):
        """Return density-free K,H in R_a=-K*na*nb*u and dEa=R_a.Vcm+H*na*nb*(Tb-Ta)."""
        ta, tb, speed = map(float, (temperature_a_ev, temperature_b_ev, drift_speed))
        if not np.isfinite([ta, tb, speed]).all() or min(ta, tb, speed) < 0:
            raise ValueError('Finite nonnegative temperatures and drift speed required')
        ma, mb = self.mass_a, self.mass_b
        total, reduced = ma+mb, ma*(mb/(ma+mb))
        variance = EV*(ta/ma+tb/mb)
        if variance == 0:
            if speed == 0:
                return dict(drag=0., exchange=0., probability_outside_table=0.)
            if not self.speed_edges[0] <= speed <= self.speed_edges[-1]:
                raise ValueError('Cold relative beam lies outside the supplied table')
            sigma = np.interp(.5*ma*speed**2/EV, self.section.energy_ev, self.section.sigma_m2)
            return dict(drag=float(reduced*speed*sigma), exchange=0., probability_outside_table=0.)
        thermal = np.sqrt(variance)
        # Resolve a narrow shifted Gaussian without removing any table interval.
        added = speed+thermal*np.arange(-12., 13.)
        added = added[(added > self.speed_edges[0]) & (added < self.speed_edges[-1])]
        edges = np.unique(np.r_[self.speed_edges, added])
        width = np.diff(edges)
        g = edges[:-1, None]+.5*width[:, None]*(1+self.nodes)
        weight = .5*width[:, None]*self.weights/thermal
        r, d = g/thermal, speed/thermal
        argument = r*d
        angular = np.divide(-np.expm1(-2*argument), argument,
                            out=np.full_like(argument, 2.), where=argument != 0)
        langevin_over_argument = np.empty_like(argument)
        small = argument < .05
        a2 = argument[small]**2
        langevin_over_argument[small] = 1/3+a2*(-1/45+a2*(2/945+a2*(-1/4725+a2*2/93555)))
        a = argument[~small]
        langevin_over_argument[~small] = (1-1/a+2*np.exp(-2*a)/(-np.expm1(-2*a)))/a
        langevin = argument*langevin_over_argument
        probability = weight*r*r/np.sqrt(2*np.pi)*np.exp(-.5*(r-d)**2)*angular
        sigma = np.interp(.5*ma*g*g/EV, self.section.energy_ev, self.section.sigma_m2)
        common = probability*sigma*g*g
        k = reduced*np.sum(common*g/variance*langevin_over_argument)
        h = reduced/(total*variance)*np.sum(common*(g-speed*langevin))
        outside = (ncx2.cdf(self.speed_edges[0]**2/variance, 3, d*d)+
                   ncx2.sf(self.speed_edges[-1]**2/variance, 3, d*d))
        if not np.isfinite([k, h, outside]).all():
            raise ValueError('Nonfinite elastic moments; no state or coefficient repair')
        return dict(drag=float(k), exchange=float(h), probability_outside_table=float(outside))

    def moments(self, density_a, density_b, temperature_a_ev, temperature_b_ev, drift_a_minus_b):
        """Paired force [N/m3] and signed heat [W/m3], including neutral recoil."""
        u = np.asarray(drift_a_minus_b, float)
        if u.shape != (3,) or not np.isfinite(u).all():
            raise ValueError('Finite three-component relative velocity required')
        na, nb = float(density_a), float(density_b)
        if not np.isfinite([na, nb]).all() or min(na, nb) <= 0:
            raise ValueError('Positive finite occupied densities required')
        result = self.coefficients(temperature_a_ev, temperature_b_ev, np.linalg.norm(u))
        k, h = na*nb*result['drag'], na*nb*result['exchange']
        power = k*float(u@u)
        exchange = h*EV*(temperature_b_ev-temperature_a_ev)
        ma, mb = self.mass_a, self.mass_b
        result.update(force_on_a=-k*u, force_on_b=k*u,
            heat_to_a=mb/(ma+mb)*power+exchange,
            heat_to_b=ma/(ma+mb)*power-exchange, drag_power=power)
        return result
