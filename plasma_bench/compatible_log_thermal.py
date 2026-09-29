"""Physical logarithmic entropy and its paired finite-mass pressure work."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E
from .entropy_scalar import EntropyScalar


@dataclass(frozen=True)
class LogThermalPair:
    entropy_rate: np.ndarray
    electric_reduced: np.ndarray
    potential: np.ndarray
    electric_full: np.ndarray
    energy_per_particle: np.ndarray
    residuals: dict


class CompatibleLogThermal:
    """Current-driven component and full density derivative; not a moving solver."""
    def __init__(self,current_model,*,gamma=5/3,density_scale=1e20,scalar=None):
        if not np.isfinite([gamma,density_scale]).all() or gamma<=1 or density_scale<=0:
            raise ValueError('Finite gamma > 1 and positive density scale required')
        self.model=current_model;self.geometry=current_model.geometry;self.mesh=current_model.mesh
        self.scalar=scalar if scalar is not None else EntropyScalar(current_model.density_map)
        if self.scalar.density_map is not current_model.density_map:
            raise ValueError('Thermal and current energy must share the same density map')
        self.gamma=gamma;self.density_scale=density_scale
        self.log_a=(gamma-1)*np.log(current_model.density/density_scale)
        self.log_a.setflags(write=False)

    def fields(self,entropy):
        s=self.scalar.sample(entropy);n=self.model.density
        with np.errstate(over='ignore',under='ignore',invalid='ignore'):
            w=E*np.exp(s)
            pressure=n*E*np.exp(self.log_a+s)
            internal=pressure/(self.gamma-1)
        if not np.isfinite([w,pressure,internal]).all() or np.any(w<=0) or np.any(internal<=0):
            raise ValueError('Thermal state exceeds positive finite numerical range; no floor')
        return dict(entropy=s,proxy_per_particle=w,pressure=pressure,internal=internal)

    def pair(self,current,entropy):
        s=self.scalar.validate(entropy);f=self.fields(s);geo=self.geometry
        current=np.asarray(current,float)
        if current.shape!=(self.model.mass.shape[0],) or not np.isfinite(current).all():
            raise ValueError('Finite compatible current coefficients required')
        j=geo.sample_electric(current);gradient=self.scalar.gradient(s)
        advective=np.sum(j*gradient,axis=-1)/E
        rate=self.scalar.solve(self.scalar.adjoint(advective))
        # dH/ds=N.T W epsilon=M_n h, hence the same physical projection must
        # enter both thermal work and pressure. It is not a fitted temperature.
        h=self.scalar.project(f['internal']/self.model.density)
        anchor=float(h.values[0,0]);small=self.scalar.sample(h.values-anchor)
        reduced=geo.project_electric(small[...,None]*gradient/E)
        potential=(anchor*s-self.gamma*h.values)/E
        full=reduced+geo.gradient(potential)
        return LogThermalPair(rate.values,reduced,potential,full,h.values,
            dict(entropy_mass=rate.residual,thermal_projection=h.residual,
                 entropy_iterations=rate.iterations,projection_iterations=h.iterations))

    def energy_density(self,entropy):return self.model.average(self.fields(entropy)['internal'])

    def energy_rate(self,entropy,entropy_rate,density_rate):
        f=self.fields(entropy);st=self.scalar.sample(entropy_rate)
        ratio=self.model.density_map.sample(density_rate)/self.model.density
        thermal=self.model.average(f['internal']*st)
        density=self.model.average(self.gamma*f['internal']*ratio)
        return dict(entropy=thermal,density=density,total=thermal+density)

    def density_energy_gradient(self,entropy):
        return self.model.density_map.adjoint(self.gamma*self.fields(entropy)['internal']/self.model.density)

    def entropy_rate(self,entropy,entropy_rate,density_rate):
        s=self.scalar.sample(entropy);st=self.scalar.sample(entropy_rate)
        nt=self.model.density_map.sample(density_rate)
        return self.model.average(self.model.density*st+nt*s)

    def proxy_rate(self,entropy,entropy_rate,density_rate):
        f=self.fields(entropy);st=self.scalar.sample(entropy_rate)
        nt=self.model.density_map.sample(density_rate)
        return self.model.average(f['proxy_per_particle']*(self.model.density*st+nt))
