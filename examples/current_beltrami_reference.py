"""Exact finite-amplitude finite-mass helical wave with constrained bulk flow."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.special import j0, j1, jn_zeros


@dataclass(frozen=True)
class CurrentBeltrami:
    radius: float = .035
    length: float = .24
    number_density: float = 1e20
    ion_mass: float = 39.948*atomic_mass
    bias: float = .005
    helicity: int = 1
    axial_mode: int = 1
    radial_mode: int = 1
    amplitude: float = 1e-6  # Vector potential, T m.
    model: str = 'two_fluid'
    electron_ev: float = 5.
    ion_ev: float = 1.

    def __post_init__(self):
        if self.helicity not in (-1, 1) or self.axial_mode == 0 or self.radial_mode < 1 or self.model not in ('two_fluid', 'hall'):
            raise ValueError('Helicity +/-1, nonzero axial mode, positive radial mode and explicit mass model required')
        if not all(np.isfinite(v) and v > 0 for v in (self.radius, self.length, self.number_density, self.ion_mass, self.electron_ev, self.ion_ev)):
            raise ValueError('Positive physical domain/species data required')

    @property
    def q(self):
        return jn_zeros(1, self.radial_mode)[-1]/self.radius

    @property
    def k(self):
        return 2*np.pi*self.axial_mode/self.length

    @property
    def total_k(self):
        return np.hypot(self.q, self.k)

    @property
    def masses(self):
        me = ME if self.model == 'two_fluid' else 0.
        return me/(self.ion_mass+me), self.ion_mass/(self.ion_mass+me), self.ion_mass*me/(self.ion_mass+me)

    @property
    def omega(self):
        alpha, beta, reduced = self.masses
        return (beta-alpha)*self.helicity*self.bias*self.k*self.total_k/(MU0*self.number_density*E*(1+reduced*self.total_k**2/(MU0*self.number_density*E*E)))

    def fields(self, r, z, time=0.):
        r, z = np.broadcast_arrays(r, z)
        phase = self.amplitude*np.exp(1j*(self.k*z-self.omega*time))
        jzero, jone = j0(self.q*r), j1(self.q*r)
        slope = np.divide(jone, r, out=np.full_like(r, self.q/2, dtype=float), where=r != 0)
        jprime = self.q*jzero-slope
        phasor = np.stack((-1j*self.k*jone*phase, self.helicity*self.total_k*jone*phase, self.q*jzero*phase), axis=-1)
        radial = np.stack((-1j*self.k*jprime*phase, self.helicity*self.total_k*jprime*phase, -self.q*self.q*jone*phase), axis=-1).real
        axial = (1j*self.k*phasor).real
        magnetic = phasor.real; bdot = (-1j*self.omega*phasor).real
        magnetic = magnetic+np.array([0., 0., self.bias])
        factor = self.helicity*self.total_k/(MU0*self.number_density*E)
        w, wr, wz, wt = factor*phasor.real, factor*radial, factor*axial, factor*bdot
        alpha, beta, reduced = self.masses
        # Cylindrical convective acceleration includes the rotating basis.
        advective = w[..., 0, None]*wr+w[..., 2, None]*wz
        advective[..., 0] -= np.divide(w[..., 1]**2, r, out=np.zeros_like(r), where=r != 0)
        advective[..., 1] += np.divide(w[..., 0]*w[..., 1], r, out=np.zeros_like(r), where=r != 0)
        electric = (beta-alpha)*np.cross(w, magnetic)+reduced/E*(wt-(beta-alpha)*advective)
        current = self.number_density*E*w
        reaction = reduced*self.number_density*advective-np.cross(current, magnetic)
        return dict(magnetic=magnetic, magnetic_perturbation=phasor.real, magnetic_rate=bdot, current=current,
                    drift=w, drift_rate=wt, drift_convective=advective, electric=electric, reaction=reaction,
                    ion_velocity=alpha*w, electron_velocity=-beta*w,
                    psi=.5*self.bias*r*r+(r*jone*phase).real)

    def initial_fields(self, mesh, time=0., order=16):
        r, z = mesh.corners; psi = self.fields(r, z, time)['psi']
        nodes, weights = np.polynomial.legendre.leggauss(order)
        g = np.zeros((mesh.nr, mesh.nz))
        for a, wa in zip(nodes, weights):
            r = (mesh.r+.5*mesh.dr*a)[:, None]
            for b, wb in zip(nodes, weights):
                field = self.fields(r, mesh.z[None]+.5*mesh.dz*b, time)['magnetic']
                g += wa*wb*mesh.dr/(4*mesh.radial_measure[:, None])*field[..., 1]
        return psi, g
