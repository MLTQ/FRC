"""Strict protocol, history, archive, continuum-flux and physical-energy audit."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from .warm_hall_cases import protocol
from .warm_hall_records import provenance
from .warm_hall_field_audit import coupled_endpoints
from .warm_audit import direct_density
from .compatible_hall_audit import magnetic_product


def isolated_endpoints(record, arrays, mesh):
    gamma = record['configuration']['reference']['gamma']
    n = arrays['number_density']; volume = arrays['dual_volumes']
    for prefix in ('initial', 'final'):
        u, g = arrays[prefix+'_entropy'], arrays[prefix+'_g']
        pe = (n/1e20)**(gamma-1)*u
        pcell = .25*(pe[:-1]+pe[1:]+np.roll(pe[:-1]+pe[1:], -1, axis=1))
        s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
        magnetic = .5*s*g*g/MU0
        ion = arrays[prefix+'_total_energy']-magnetic-pcell/(gamma-1)
        assert np.all(ion > 0) and np.all(pe > 0)
        values = dict(magnetic_energy_j=np.sum(mesh.volumes*magnetic), electron_internal_energy_j=np.sum(volume*pe)/(gamma-1),
                      ion_internal_energy_j=np.sum(mesh.volumes*ion), total_energy_j=np.sum(mesh.volumes*arrays[prefix+'_total_energy']),
                      entropy_proxy_j=np.sum(volume*u), entropy_log_particles=np.sum(volume*n*np.log(u/(n*E))))
        for k, v in values.items():
            np.testing.assert_allclose(v, record[prefix][k], rtol=5e-13, atol=1e-25)
        if prefix == 'final':
            scale = (gamma-1)*(abs(arrays[prefix+'_total_energy'])+abs(magnetic))+abs(pcell)
            assert np.all(abs((gamma-1)*ion-arrays['final_ion_pressure']) <= 32*np.finfo(float).eps*scale)


def history(record, arrays):
    h, cfg = record['history_columns'], record['configuration']
    assert all(len(v) == record['accepted_steps']+1 and np.isfinite(v).all() for v in h.values())
    for k, v in h.items():
        assert v[0] == record['initial'][k] and v[-1] == record['final'][k]
    time = np.asarray(h['time_s']); end = record['end_time_s']
    assert time[0] == 0 and np.all(np.diff(time) > 0) and abs(time[-1]-end) <= 2e-13*end
    if cfg['steps'] is not None:
        assert record['accepted_steps'] == cfg['steps']
        np.testing.assert_allclose(np.diff(time), end/cfg['steps'], rtol=5e-13)
    for key in ('min_ion_pressure_pa', 'min_electron_pressure_pa'):
        assert min(h[key]) > 0
    coupled = cfg['kind'] in ('wave', 'nonlinear')
    if coupled:
        assert h['steps'] == list(range(record['accepted_steps']+1)) and time[-1] == record['actual_time_s']
        assert not np.any(h['physical_ohmic_heat_j'])
        assert record['final']['max_hall_nonlinear_residual'] <= cfg['tolerance']
        assert record['final']['max_pressure_nonlinear_residual'] <= cfg['tolerance']
        cap = end*cfg['phase_step']/cfg['phase'] if cfg['kind'] == 'wave' else cfg['dt_max']
        assert max(np.diff(time)) <= cap*(1+5e-13)
    else:
        assert max(h['nonlinear_residual']) <= cfg['tolerance']
    thermal = np.array(h['electron_internal_energy_j'])+h['ion_internal_energy_j']
    internal = h['internal_energy_j'] if coupled else np.array(h['total_energy_j'])-h['magnetic_energy_j']
    np.testing.assert_allclose(thermal, internal, rtol=5e-15)
    keys = ('mass_kg', 'total_energy_j', 'electron_entropy_proxy_integral_j') if coupled else ('total_energy_j', 'entropy_proxy_j')
    result = {f'max_{k}_relative_drift': float(np.max(abs(np.array(h[k])-h[k][0]))/abs(h[k][0])) for k in keys}
    sk = 'electron_log_entropy_particles' if coupled else 'entropy_log_particles'
    particles = (record['initial']['mass_kg']/(cfg['reference']['ion_mass'] if cfg['kind'] == 'wave' else cfg['initial']['ion_mass'])) if coupled else np.sum(arrays['dual_volumes']*arrays['number_density'])
    changes = (np.array(h[sk])-h[sk][0])/particles
    result.update(min_log_entropy_change_per_particle=float(changes.min()), max_log_entropy_change_per_particle=float(changes.max()))
    return result


def audit(directory='outputs/warm_hall_verification', configurations=None, *, assembled=True):
    directory = Path(directory)
    configurations = protocol() if configurations is None else configurations
    data = json.loads((directory/'measurements.json').read_text()) if assembled else dict(
        completed=True, provenance=provenance(), cases=[json.loads((directory/(c['name']+'.json')).read_text()) for c in configurations])
    assert data['completed'] and data['provenance'] == provenance()
    assert [r['configuration'] for r in data['cases']] == configurations
    records, arrays, meshes, results = {}, {}, {}, {}
    for record in data['cases']:
        cfg = record['configuration']; name = cfg['name']; ref = cfg['reference']
        assert record['completed'] and record['provenance'] == data['provenance'] and name not in records
        assert record == json.loads((directory/(name+'.json')).read_text())
        path = directory/record['fields']['file']
        assert path.name == name+'.npz' and hashlib.sha256(path.read_bytes()).hexdigest() == record['fields']['sha256']
        with np.load(path, allow_pickle=False) as opened:
            a = {k: opened[k] for k in opened.files}
        assert all(np.isfinite(v).all() for v in a.values())
        mesh = MHDMesh(ref['radius'], ref['length']/2, cfg['nr'], cfg['nz'])
        for key, expected in (('radial_faces', mesh.rf), ('axial_faces', mesh.zf), ('primal_volumes', mesh.volumes)):
            np.testing.assert_array_equal(a[key], expected)
        if 'dual_volumes' in a:
            _, volume = direct_density(mesh, np.ones((mesh.nr, mesh.nz, 1)), 1.)
            np.testing.assert_allclose(a['dual_volumes'], volume, rtol=8*mesh.nr*np.finfo(float).eps)
        result = {}
        if cfg['kind'] == 'pressure_spatial':
            rate = -(np.diff(a['radial_flux'], axis=0)+np.diff(a['axial_flux'], axis=1))/mesh.volumes
            np.testing.assert_array_equal(rate, a['energy_rate'])
            exact_div = -(np.diff(a['reference_radial_flux'], axis=0)+np.diff(a['reference_axial_flux'], axis=1))/mesh.volumes
            scale = np.linalg.norm(a['reference_pressure_work'])
            result['independent_continuum_energy_flux_relative_residual'] = float(np.linalg.norm(exact_div-a['reference_energy_rate'])/scale)
            assert result['independent_continuum_energy_flux_relative_residual'] < 2e-11
        elif cfg['kind'] == 'trace_diffusion':
            h = []
            for prefix in ('initial', 'final'):
                psi, g = a[prefix+'_psi'], a[prefix+'_g']
                h.append(.5*np.sum(mesh.volumes*magnetic_product(mesh, psi, g, psi, g)))
            error = (h[1]-h[0]+record['physical_heat_j'])/h[0]
            result['independent_energy_heat_relative_residual'] = float(error)
            assert abs(error) < 2e-12 and record['physical_heat_j'] > 0
        else:
            if cfg['kind'] == 'pressure_evolve':
                isolated_endpoints(record, a, mesh)
            else:
                result.update(coupled_endpoints(record, a, mesh))
            result.update(history(record, a))
        records[name], arrays[name], meshes[name], results[name] = record, a, mesh, result
    return data, records, arrays, meshes, results
