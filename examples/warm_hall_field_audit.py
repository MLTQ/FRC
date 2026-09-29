"""Independent physical warm-wave pressure and energy replay at higher quadrature."""
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.toroidal_state import ToroidalState, primitive
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction
from .warm_audit import direct_density
from .compatible_hall_audit import magnetic_product
from .warm_hall_norms import sample_q1
from .warm_hall_reference import WarmHallWave
from .compatible_hall_references import wave_state


def replay(mesh, ref, arrays, prefix, diagnostic):
    u, psi, entropy = (arrays[prefix+'_'+k] for k in ('conserved', 'psi', 'entropy'))
    n, volume = direct_density(mesh, u, ref.ion_mass)
    # Subtracting squared neighboring radii loses O(nr*eps), unlike nr=4 controls.
    np.testing.assert_allclose(volume, arrays['dual_volumes'], rtol=8*mesh.nr*np.finfo(float).eps)
    pe = (n/1e20)**(ref.gamma-1)*entropy
    s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
    kinetic = .5*(u[..., 1]**2+u[..., 2]**2+u[..., 4]**2/s)/u[..., 0]
    magnetic = .5*magnetic_product(mesh, psi, u[..., 5], psi, u[..., 5])
    internal = u[..., 3]-kinetic-magnetic
    pcell = .25*(pe[:-1]+pe[1:]+np.roll(pe[:-1]+pe[1:], -1, axis=1))
    pi = (ref.gamma-1)*internal-pcell
    assert np.all(n > 0) and np.all(pe > 0) and np.all(pi > 0)
    fields = dict(mass_kg=u[..., 0], total_energy_j=u[..., 3], kinetic_energy_j=kinetic,
                  magnetic_energy_j=magnetic, internal_energy_j=internal, ion_internal_energy_j=pi/(ref.gamma-1))
    values = {k: np.sum(mesh.volumes*v) for k, v in fields.items()}
    values.update(electron_internal_energy_j=np.sum(volume*pe)/(ref.gamma-1),
                  electron_entropy_proxy_integral_j=np.sum(volume*entropy),
                  electron_log_entropy_particles=np.sum(volume*n*np.log(entropy/(n*E))))
    for key, value in values.items():
        np.testing.assert_allclose(value, diagnostic[key], rtol=4e-12, atol=1e-28, err_msg=key)
    return pe, pi


def pressure_norms(mesh, archive, ref, time, order):
    """Integrate physical errors in x,z; do not use electric-space samples."""
    state = ToroidalState(archive['final_conserved'], archive['final_psi'], 'compatible')
    rec = ToroidalReconstruction(mesh, primitive(mesh, state, ref.gamma), order=2)
    n, _ = direct_density(mesh, state.conserved, ref.ion_mass)
    pe = (n/1e20)**(ref.gamma-1)*archive['final_entropy']
    nodes, weights = np.polynomial.legendre.leggauss(order)
    lo, hi = .5*mesh.rf[:-1]**2, .5*mesh.rf[1:]**2
    numerator = np.zeros(4); denominator = np.zeros(4)
    minimum_pi = np.inf
    for a, wa in zip(nodes, weights):
        x = .5*(lo+hi+a*(hi-lo))[:, None]
        r = np.sqrt(2*x)
        for b, wb in zip(nodes, weights):
            z = mesh.z[None]+.5*mesh.dz*b
            q = rec.primitive(rec.raw(r[:, 0], .5*b), rec.auxiliary(r[:, 0], .5*b), r)
            actual_pe, actual_n = sample_q1(mesh, pe, x, z), sample_q1(mesh, n, x, z)
            actual_pi = q[..., 3]-actual_pe
            minimum_pi = min(minimum_pi, float(actual_pi.min()))
            exact = ref.fields(r, z, time)
            true_pe = ref.electron_pressure+exact['electron_pressure']
            true_pi = ref.pressure-ref.electron_pressure+exact['ion_pressure']
            true_n = ref.number_density+exact['density']/ref.ion_mass
            values = (actual_pe, actual_pi, actual_pe/(actual_n*E), actual_pi/(q[..., 0]/ref.ion_mass*E))
            target = (true_pe, true_pi, true_pe/(true_n*E), true_pi/(true_n*E))
            variation = (exact['electron_pressure'], exact['ion_pressure'],
                         target[2]-ref.electron_pressure/(ref.number_density*E),
                         target[3]-(ref.pressure-ref.electron_pressure)/(ref.number_density*E))
            for k in range(4):
                w = .25*wa*wb*mesh.volumes
                numerator[k] += np.sum(w*(values[k]-target[k])**2)
                denominator[k] += np.sum(w*variation[k]**2)
    keys = ('electron_pressure', 'ion_pressure', 'electron_temperature', 'ion_temperature')
    return dict(relative_errors=dict(zip(keys, np.sqrt(numerator/denominator).tolist())),
                rms_errors=dict(zip(keys, np.sqrt(numerator/np.sum(mesh.volumes)).tolist())),
                minimum_reconstructed_ion_pressure_pa=minimum_pi)


def coupled_endpoints(record, arrays, mesh):
    from types import SimpleNamespace
    case = record['configuration']; wave = case['kind'] == 'wave'
    ref = WarmHallWave(**case['reference']) if wave else SimpleNamespace(
        ion_mass=case['initial']['ion_mass'], gamma=case['reference']['gamma'])
    for prefix in ('initial', 'final'):
        pe, pi = replay(mesh, ref, arrays, prefix, record[prefix])
        np.testing.assert_allclose(pe, arrays[prefix+'_electron_pressure'], rtol=5e-15)
        if prefix == 'final':
            np.testing.assert_allclose(pi, arrays['final_ion_pressure'], rtol=5e-13)
        psi = arrays[prefix+'_psi']
        np.testing.assert_array_equal(psi[:, 0], psi[:, -1])
        np.testing.assert_allclose(psi[0], 0., atol=1e-20, rtol=0)
        np.testing.assert_allclose(psi[-1], psi[-1, 0], atol=1e-18, rtol=0)
    result = {}
    if wave:
        exact = wave_state(mesh, ref, order=18)
        np.testing.assert_allclose(arrays['initial_conserved'], exact.conserved, rtol=2e-13, atol=1e-24)
        np.testing.assert_array_equal(arrays['initial_psi'], exact.psi)
        norm8 = pressure_norms(mesh, arrays, ref, record['actual_time_s'], 8)
        norm12 = pressure_norms(mesh, arrays, ref, record['actual_time_s'], 12)
        result['independent_physical_norms'] = norm12
        result['quadrature_8_to_12_max_relative_norm_change'] = max(abs(norm8['relative_errors'][k]-v) for k, v in norm12['relative_errors'].items())
        result['measurement_to_independent_max_relative_norm_change'] = max(abs(record['metrics'][k+'_relative_perturbation_L2']-v) for k, v in norm12['relative_errors'].items())
        assert norm12['minimum_reconstructed_ion_pressure_pa'] > 0
    return result
