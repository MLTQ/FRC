"""Explicit lumped graph pressure/work bridge, with no high-order accuracy claim."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E
from .compatible_charge import CompatibleCharge
from .electron_pressure import ElectronPressure
from .thermal_means import entropy_mean


@dataclass(frozen=True)
class PressurePair:
    entropy_rate: np.ndarray
    electric_reduced: np.ndarray
    potential: np.ndarray
    electric_full: np.ndarray


class CompatibleSpeciesPressure:
    """At supplied positive nodal n, pressure/current thermal exchange only."""
    def __init__(self, geometry, nodal_density, *, gamma=5/3, density_scale=1e20):
        self.field_geometry=geometry;self.geometry=CompatibleCharge(geometry)
        # Reuse precisely the existing lumped thermodynamic definition, not its
        # old magnetic mass or charge stream.
        self.thermo=ElectronPressure(geometry.mesh,nodal_density,gamma=gamma,
                                    density_scale=density_scale)

    def pair(self, g, entropy_density):
        geo=self.geometry;thermo=self.thermo
        u=geo.scalar(entropy_density,'species entropy density')
        if np.any(u<=0):raise ValueError('Positive entropy density required; no repair')
        w=u/thermo.density
        wr=entropy_mean(w[1:],w[:-1]);wz=entropy_mean(w,np.roll(w,-1,axis=1))
        ir,iz=geo.current_flux(g)
        flux=np.zeros_like(ir);flux[1:-1]=ir[1:-1]*wr/E
        rate=geo.divergence(flux,iz*wz/E)
        anchor=np.sum(geo.volumes*u)/np.sum(geo.volumes*thermo.density)
        ar,az=geo.gradient(thermo.energy_weight)
        reduced=np.r_[(-(wr-anchor)*ar/E).ravel(),
                      np.zeros((geo.mesh.nr-1)*geo.mesh.nz),
                      (-(wz-anchor)*az/E).ravel()]
        potential=-(thermo.a*w+anchor*thermo.energy_weight)/E
        full=reduced+self.field_geometry.gradient(potential)
        return PressurePair(rate,reduced,potential,full)

    def thermal_energy_rate(self, entropy_rate):
        return float(np.sum(self.geometry.volumes*self.thermo.energy_weight*
                            self.geometry.scalar(entropy_rate)))
