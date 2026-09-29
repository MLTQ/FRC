"""Finite-species-mass current/thermal subsystem at prescribed n and zero bulk flow."""
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME
from scipy.sparse import block_diag
from .hall_trace import TraceHallCurl
from .hall_matrices import magnetic_hessian, transport_matrix
from .electron_pressure import ElectronPressure
from .dual_mass import DualMass
from .thermal_means import entropy_mean


class CurrentThermal:
    """Zero barycentric flow constraint; not a complete moving two-fluid solver."""
    def __init__(self, mesh, number_density, ion_mass, *, model='two_fluid', gamma_e=5/3, gamma_i=5/3,
                 density_scale=1e20, quadrature_order=5):
        if not np.isfinite(ion_mass) or ion_mass <= 0 or model not in ('two_fluid', 'hall'):
            raise ValueError('Positive physical ion mass and explicit two_fluid or hall model required')
        self.mesh, self.geometry = mesh, TraceHallCurl(mesh, quadrature_order)
        self.electron = ElectronPressure(mesh, number_density, gamma=gamma_e, density_scale=density_scale)
        self.ion = ElectronPressure(mesh, number_density, gamma=gamma_i, density_scale=density_scale)
        self.mapping = DualMass(mesh)
        self.ion_mass, self.electron_mass = ion_mass, ME if model == 'two_fluid' else 0.
        total_mass = self.ion_mass+self.electron_mass
        self.alpha, self.beta = self.electron_mass/total_mass, self.ion_mass/total_mass
        self.delta, self.reduced_mass = self.beta-self.alpha, self.ion_mass*self.alpha
        self.density = self.geometry.scalar_samples(self.electron.density)
        self.inertia_coefficient = self.reduced_mass/(self.density*E*E)
        self.density.setflags(write=False)
        self.inertia_coefficient.setflags(write=False)
        self.mass = block_diag(self.geometry.mass, format='csc')
        self.hessian = magnetic_hessian(mesh)
        self.inertia = transport_matrix(self.geometry, self.geometry.magnetic_samples(
            np.zeros((mesh.nr+1, mesh.nz+1)), np.zeros((mesh.nr, mesh.nz))), None, self.inertia_coefficient)

    def current(self, packed):
        source = self.geometry.split_electric(-self.geometry.faraday.T @ (self.hessian @ packed))
        return np.concatenate([solve(v) for solve, v in zip(self.geometry.solves, source)])

    def star_field(self, packed, current):
        """B*=B+curl(mu J/(n e²)), including convective current inertia."""
        drift = self.geometry.sample_electric(current)/(self.density[..., None]*E)
        coefficients = self.geometry.project_electric(drift)
        star = packed-self.reduced_mass/E*(self.geometry.faraday @ coefficients)
        return self.geometry.magnetic_samples(*self.geometry.unpack(star))

    def gradient(self, scalar):
        er, ez = self.electron.geometry.gradient(scalar)
        return np.concatenate((er.ravel(), np.zeros(self.geometry.maps[1].shape[1]), ez.ravel()))

    def reduced_pressure(self, species, entropy):
        """Separate a complete discrete gradient before the stiff curl solve."""
        geo = species.geometry; u = geo.scalar(entropy)
        if np.any(u <= 0):
            raise ValueError('Positive species entropy required; no floor applied')
        w = u/species.density
        anchor = np.sum(geo.volumes*u)/np.sum(geo.volumes*species.density)
        wr, wz = entropy_mean(w[1:], w[:-1]), entropy_mean(w, np.roll(w, -1, axis=1))
        ar, az = geo.gradient(species.energy_weight)
        er, ez = -(wr-anchor)*ar/E, -(wz-anchor)*az/E
        electric = np.concatenate((er.ravel(), np.zeros(self.geometry.maps[1].shape[1]), ez.ravel()))
        potential = -(species.a*w+anchor*species.energy_weight)/E
        return electric, potential

    def pressure(self, ue, ui):
        ee, phie = self.reduced_pressure(self.electron, ue)
        ei, phii = self.reduced_pressure(self.ion, ui)
        return self.beta*ee-self.alpha*ei, self.beta*phie-self.alpha*phii

    def thermal_rates(self, g, ue, ui):
        return self.beta*self.electron.rhs(g, ue)[1], -self.alpha*self.ion.rhs(g, ui)[1]

    def energy_density(self, packed, ue, ui):
        psi, g = self.geometry.unpack(packed)
        current = self.geometry.sample_electric(self.current(packed))
        kinetic = np.sum(self.geometry.weight*self.inertia_coefficient*np.sum(current*current, axis=-1), axis=(2, 3))/(2*self.mesh.volumes)
        electron = self.mapping.primal_average(self.electron.energy_weight*ue)
        ion = self.mapping.primal_average(self.ion.energy_weight*ui)
        return dict(magnetic=self.geometry.energy_density(psi, g), current_kinetic=kinetic, electron=electron, ion=ion)

    def diagnostics(self, packed, ue, ui):
        densities = self.energy_density(packed, ue, ui)
        out = {k+'_energy_j': float(np.sum(self.mesh.volumes*v)) for k, v in densities.items()}
        out['total_energy_j'] = sum(out.values())
        for label, species, u in (('electron', self.electron, ue), ('ion', self.ion, ui)):
            d = species.diagnostics(self.geometry.unpack(packed)[1], u)
            out[label+'_log_entropy_particles'] = d['entropy_log_particle_units']
            out[label+'_entropy_proxy_j'] = d['entropy_proxy_integral_j']
            out['min_'+label+'_pressure_pa'] = d['min_electron_pressure_pa']
        return out
