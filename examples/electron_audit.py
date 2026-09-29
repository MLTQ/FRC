"""Audit electron records and independently replay physical energy and entropy."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from .electron_cases import protocol
from .electron_records import provenance, relative


def load(directory=Path('outputs/electron_pressure_verification')):
    collection = json.loads((directory/'measurements.json').read_text())
    assert collection['completed'] and collection['provenance'] == provenance()
    assert [r['configuration'] for r in collection['cases']] == protocol()
    records, archives, audits = {}, {}, {}
    for record in collection['cases']:
        case = record['configuration']; name = case['name']
        assert name not in records and record['completed'] and record['provenance'] == collection['provenance']
        assert record == json.loads((directory/(name+'.json')).read_text())
        path = directory/record['fields']['file']
        assert path.name == name+'.npz' and hashlib.sha256(path.read_bytes()).hexdigest() == record['fields']['sha256']
        with np.load(path, allow_pickle=False) as data:
            a = {k: data[k] for k in data.files}
        assert all(np.isfinite(v).all() for v in a.values())
        nr, nz = case['nr'], case['nz']; ref = case['reference']
        rf, zf = np.linspace(0, ref['radius'], nr+1), np.linspace(-ref['length']/2, ref['length']/2, nz+1)
        np.testing.assert_array_equal(a['radial_faces'], rf); np.testing.assert_array_equal(a['axial_faces'], zf)
        x, dz = rf*rf/2, ref['length']/nz
        mass = a['primal_volumes']*(x[:-1]+x[1:])[:, None]/MU0
        dual_faces = np.r_[x[0], (x[:-1]+x[1:])/2, x[-1]]
        np.testing.assert_allclose(a['dual_volumes'], np.broadcast_to(2*np.pi*np.diff(dual_faces)[:, None]*dz, (nr+1, nz)), rtol=2e-15)
        np.testing.assert_allclose(a['primal_volumes'], np.broadcast_to(2*np.pi*np.diff(x)[:, None]*dz, (nr, nz)), rtol=2e-15)
        assert a['number_density'].shape == (nr+1, nz) and np.all(a['number_density'] > 0)
        if record['kind'] == 'spatial':
            ir, iz = a['current_radial'], a['current_axial']
            residual = np.diff(ir, axis=0)+iz-np.roll(iz, 1, axis=1)
            charge = float(abs(residual).max()/max(abs(ir).max(), abs(iz).max()))
            thermal = a['dual_volumes']*(a['number_density']/case['density_scale'])**(ref['gamma']-1)/(ref['gamma']-1)*a['thermal_rate']
            magnetic = mass*a['g']*a['battery_rate']
            work = float(abs(thermal.sum()+magnetic.sum())/(np.sum(abs(thermal))+np.sum(abs(magnetic))))
            entropy = a['dual_volumes']*a['number_density']/a['entropy_density']*a['thermal_rate']
            entropy_error = float(abs(entropy.sum())/np.sum(abs(entropy)))
            audits[name] = dict(charge_backward_error=charge, work_backward_error=work, logarithmic_entropy_rate_backward_error=entropy_error)
            assert max(charge, work, entropy_error) < 2e-13
        else:
            audits[name] = trajectory(record, a, mass)
        records[name], archives[name] = record, a
    return records, archives, audits


def trajectory(record, a, magnetic_mass):
    case = record['configuration']; gamma = case['reference']['gamma']
    volume, n = a['dual_volumes'], a['number_density']
    particles = float(np.sum(volume*n))
    for prefix in ('initial', 'final'):
        p, u, g = a[prefix+'_pressure'], a[prefix+'_entropy_density'], a[prefix+'_g']
        assert p.shape == u.shape == volume.shape and g.shape == magnetic_mass.shape
        assert np.all(p > 0) and np.all(u > 0)
        np.testing.assert_allclose(p, (n/case['density_scale'])**(gamma-1)*u, rtol=3e-15)
        magnetic = float(.5*np.sum(magnetic_mass*g*g))
        thermal = float(np.sum(volume*p/(gamma-1)))
        entropy = float(np.sum(volume*n*(np.log(p/(n*E))-(gamma-1)*np.log(n/case['density_scale']))))
        d = record[prefix]
        assert abs(magnetic-d['magnetic_energy_j']) <= 5e-15*max(magnetic, 1e-300)
        assert abs(thermal-d['electron_internal_energy_j']) < 5e-15*thermal
        assert abs(entropy-d['entropy_log_particle_units'])/particles < 2e-14
    history = record['history_columns']
    steps = record['accepted_steps']
    assert all(len(v) == steps+1 and np.isfinite(v).all() for v in history.values())
    assert history['time_s'][0] == 0 and np.all(np.diff(history['time_s']) > 0)
    assert history['time_s'][-1] == record['actual_time_s']
    assert (record['failure'] is not None) == case['expected_failure']
    if record['failure']:
        assert steps == 0 and record['actual_time_s'] == 0
        np.testing.assert_array_equal(a['initial_g'], a['final_g'])
        np.testing.assert_array_equal(a['initial_entropy_density'], a['final_entropy_density'])
    else:
        assert steps == case['steps'] and record['actual_time_s'] == record['end_time_s']
    assert min(history['min_electron_pressure_pa']) > 0 and record['max_nonlinear_residual'] <= case['tolerance']
    energy = np.asarray(history['total_energy_j'])
    entropy = np.asarray(history['entropy_log_particle_units'])
    proxy = np.asarray(history['entropy_proxy_integral_j'])
    return dict(max_energy_relative_drift=float(max(abs(energy-energy[0]))/energy[0]),
                max_logarithmic_entropy_change_per_particle=float(max(abs(entropy-entropy[0]))/particles),
                max_proxy_relative_drift=float(max(abs(proxy-proxy[0]))/proxy[0]), particles=particles,
                completed_time=not bool(record['failure']))


def difference(a, b, g_reference, p_reference):
    s = .5*(a['radial_faces'][:-1]**2+a['radial_faces'][1:]**2)[:, None]
    return dict(magnetic_relative_to_reference_change=relative(a['final_g']-g_reference, g_reference-b['initial_g'], a['primal_volumes']*s),
                pressure_relative_to_reference_change=relative(a['final_pressure']-p_reference, p_reference-b['initial_pressure'], a['dual_volumes']))


def comparisons(records, arrays):
    times, meshes, controls = {}, {}, {}
    for family, counts, fine in [('drift', (64, 128, 256), 1024), ('nonlinear', (10, 20, 40), 320)]:
        ref = arrays[f'{family}_time_{fine}']
        for count in counts:
            name = f'{family}_time_{count}'
            times[name] = difference(arrays[name], ref, ref['final_g'], ref['final_pressure'])
    for family in ('battery', 'nonlinear'):
        ref = arrays[f'{family}_n64']
        for nr in (8, 16, 32):
            name = f'{family}_n{nr}'; a = arrays[name]; ratio = 64//nr
            def average(g):
                return (g*ref['primal_volumes']).reshape(nr, ratio, 2*nr, ratio).sum(axis=(1, 3))/a['primal_volumes']
            g, p = average(ref['final_g']), ref['final_pressure'][::ratio, ::ratio]
            initial = dict(initial_g=average(ref['initial_g']), initial_pressure=ref['initial_pressure'][::ratio, ::ratio])
            meshes[name] = difference(a, initial, g, p)
    ref = arrays['nonlinear_n16']
    for name in ('nonlinear_scale', 'nonlinear_tolerance_1e-08', 'nonlinear_tolerance_1e-12'):
        controls[name] = difference(arrays[name], ref, ref['final_g'], ref['final_pressure'])
    return dict(time=times, nested_mesh=meshes, controls=controls,
                norm='Difference divided by reference change; g annular magnetic norm, pressure lumped dual nodal norm. Mesh g is volume-restricted, p sampled at coincident nodes; not continuum field truth.')


def physical_drift(record, archive, order=8):
    """Independent continuous Q1 pressure error between scalar nodes."""
    case = record['configuration']; ref = case['reference']
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes+1)/2
    rf, zf, p = archive['radial_faces'], archive['axial_faces'], archive['final_pressure']
    volume = archive['primal_volumes']
    velocity = -2*case['slope']/(MU0*ref['density']*E)
    error, norm = 0., 0.
    for i, a in enumerate(t):
        down = (1-a)*p[:-1]+a*p[1:]
        up = np.roll(down, -1, axis=1)
        for j, b in enumerate(t):
            z = (1-b)*zf[:-1]+b*zf[1:]
            exact = 80*(1+case['amplitude']*np.cos(2*np.pi*(z-velocity*record['actual_time_s'])/ref['length']))[None]
            actual = (1-b)*down+b*up
            weight = .25*weights[i]*weights[j]*volume
            error += np.sum(weight*(actual-exact)**2)
            norm += np.sum(weight*(exact-80)**2)
    return float(np.sqrt(error/norm))
