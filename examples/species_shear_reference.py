"""Exact entropy translation along a prescribed sheared axial current."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .regular_coupling_reference import initialize as base


def density(ref, r):
    return ref.number_density*(1+.3*(r/ref.radius)**2)


def current_z(ref, r, sign=1):
    return sign*2*ref.toroidal/MU0*(1-2*(r/ref.radius)**2)


def entropy(ref, r, z, time, fraction, sign=1):
    speed = fraction*current_z(ref, r, sign)/(E*density(ref, r))
    return np.log(2.)+.08*(r/ref.radius)**2+.12*np.cos(np.pi/ref.half_length*(z-speed*time))


def initialize(spec):
    ref, d, u, m, _ = base(spec)
    n = np.broadcast_to(density(ref, d.r), u[..., 0].shape)
    u[...] = 0.
    u[..., 0] = n
    for k, a in enumerate((-m.beta, m.alpha)):
        u[..., k+4] = n*entropy(ref, d.r, d.z, 0., a, spec.get('sign', 1))
    m = RegularCurrent(m.mapping, n, ref.ion_mass)
    mesh = m.mesh
    xmean = .25*(mesh.rf[:-1]**2+mesh.rf[1:]**2)
    g = np.broadcast_to((spec.get('sign', 1)*ref.toroidal*(1-2*xmean/ref.radius**2))[:, None],
                        (mesh.nr, mesh.nz)).copy()
    field = HarmonicMagnetic(mesh, np.zeros((mesh.nr+1, mesh.nz+1)), g, ref.axial_bias)
    return ref, d, u, m, field
