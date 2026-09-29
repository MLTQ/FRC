"""Continuous-volume species pressure and temperature error integration."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.toroidal_state import primitive
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction


def physical_norms(mesh, state, electron_pressure, dual_density, ref, time, order=8, reconstruction_order=2):
    """Compare Q1 electrons and reconstructed total-minus-electron ions in SI."""
    q = primitive(mesh, state, ref.gamma)
    reconstruction = ToroidalReconstruction(mesh, q, order=reconstruction_order)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes+1)/2
    errors = np.zeros(4); denominators = np.zeros(4)
    minimum = np.full(4, np.inf)
    baselines = np.array([ref.pe0, ref.p0-ref.pe0, ref.electron_ev, ref.ion_ev])
    for i, tx in enumerate(t):
        radius = np.sqrt((1-tx)*mesh.rf[:-1]**2+tx*mesh.rf[1:]**2)
        for j, tz in enumerate(t):
            qq = reconstruction.primitive(reconstruction.raw(radius, tz-.5), reconstruction.auxiliary(radius, tz-.5), radius[:, None])
            def sample(v):
                radial = (1-tx)*v[:-1]+tx*v[1:]
                return (1-tz)*radial+tz*np.roll(radial, -1, axis=1)
            pe, n = sample(electron_pressure), sample(dual_density)
            pi = qq[..., 3]-pe
            actual = np.stack((pe, pi, pe/(n*E), pi/(qq[..., 0]/ref.ion_mass*E)))
            exact = ref.fields(mesh.zf[:-1]+tz*mesh.dz, time)
            exact_n = exact['density']/ref.ion_mass
            target = np.stack((exact['electron_pressure'], exact['ion_pressure'],
                               exact['electron_pressure']/(exact_n*E), exact['ion_pressure']/(exact_n*E)))[:, None, :]
            measure = .25*weights[i]*weights[j]*mesh.volumes
            errors += np.sum(measure*(actual-target)**2, axis=(1, 2))
            denominators += np.sum(measure*(target-baselines[:, None, None])**2, axis=(1, 2))
            minimum = np.minimum(minimum, np.min(actual, axis=(1, 2)))
    names = ('electron_pressure', 'reconstructed_ion_pressure', 'electron_temperature', 'reconstructed_ion_temperature')
    return {name: dict(rms_error=float(np.sqrt(errors[i]/np.sum(mesh.volumes))),
                       perturbation_relative_L2=float(np.sqrt(errors[i]/denominators[i])) if denominators[i] else None,
                       min_reconstructed_value=float(minimum[i]), unit='Pa' if i < 2 else 'eV') for i, name in enumerate(names)}
