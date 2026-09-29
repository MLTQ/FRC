"""Strict protocol/archive audit and independent cylindrical energy replay."""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from .compatible_hall_cases import protocol
from .compatible_hall_controls import protocol as controls
from .compatible_hall_records import provenance


def magnetic_product(mesh, psi, g, other_psi, other_g):
    """Independent direct physical Gauss integration of B dot other_B / mu0."""
    nodes, weights = np.polynomial.legendre.leggauss(24)
    t = (nodes+1)/2
    result = np.zeros((mesh.nr, mesh.nz))
    lo, hi = mesh.rf[:-1, None], mesh.rf[1:, None]
    radial_measure = (hi*hi-lo*lo)/2
    for i in range(len(t)):
        r = lo+(hi-lo)*t[i]
        x = (r*r-lo*lo)/(hi*hi-lo*lo)
        for j in range(2):
            z = .5+(-1 if j == 0 else 1)/(2*np.sqrt(3))
            def fields(p, tor):
                br = -((1-x)*(p[:-1, 1:]-p[:-1, :-1])+x*(p[1:, 1:]-p[1:, :-1]))/(r*mesh.dz)
                bz = ((1-z)*(p[1:, :-1]-p[:-1, :-1])+z*(p[1:, 1:]-p[:-1, 1:]))/radial_measure
                return br, r*tor, bz
            a, b = fields(psi, g), fields(other_psi, other_g)
            result += weights[i]*.25*r*(hi-lo)/radial_measure*sum(v*w for v, w in zip(a, b))/MU0
    return result


def collection(directory, configurations):
    directory = Path(directory)
    data = json.loads((directory/'measurements.json').read_text())
    assert data['completed'] and data['provenance'] == provenance()
    assert [r['configuration'] for r in data['cases']] == configurations
    records, archives, meshes, audits = {}, {}, {}, {}
    for expected, record in zip(configurations, data['cases']):
        name = expected['name']
        assert name not in records and record['completed']
        assert record == json.loads((directory/(name+'.json')).read_text())
        assert record['provenance'] == data['provenance']
        path = directory/record['fields']['file']
        assert record['fields']['file'] == name+'.npz'
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['fields']['sha256']
        with np.load(path, allow_pickle=False) as opened:
            archive = {key: opened[key] for key in opened.files}
        assert all(np.isfinite(a).all() for a in archive.values())
        ref = expected.get('reference', expected)
        mesh = MHDMesh(ref['radius'], ref['length']/2, expected['nr'], expected['nz'])
        np.testing.assert_array_equal(archive['radial_faces'], mesh.rf)
        np.testing.assert_array_equal(archive['axial_faces'], mesh.zf)
        if record['kind'] == 'spatial':
            rate = magnetic_product(mesh, archive['psi'], archive['g'], archive['psi_rate'], archive['g_rate'])
            scale = record['metrics']['work_scale_w']
            rate_norm = np.linalg.norm(rate)
            replay_error = np.linalg.norm(rate-archive['magnetic_rate'])
            if rate_norm == 0:
                assert replay_error == 0
            audits[name] = dict(replayed_work_relative=abs(float(np.sum(mesh.volumes*rate)))/scale,
                                local_rate_replay_relative=float(replay_error/rate_norm) if rate_norm else 0.)
            assert audits[name]['replayed_work_relative'] < 2e-12
            assert audits[name]['local_rate_replay_relative'] < 2e-12
        else:
            audits[name] = trajectory(record, archive, mesh, ref.get('gamma', 5/3))
        records[name], archives[name], meshes[name] = record, archive, mesh
    return records, archives, meshes, audits


def trajectory(record, archive, mesh, gamma):
    history, end, final = record['history_columns'], record['end_time_s'], record['final']
    np.testing.assert_array_equal(archive['volumes'], mesh.volumes)
    assert all(len(v) == final['steps']+1 and np.isfinite(v).all() for v in history.values())
    time = np.asarray(history['time_s'])
    assert time[0] == 0 and np.all(np.diff(time) > 0)
    assert abs(time[-1]-end) <= 4*np.finfo(float).eps*end and final['time_s'] == end
    assert min(history['min_pressure_pa']) > 0
    heat = np.asarray(history['physical_ohmic_heat_j'])
    assert heat[0] == 0 and np.all(np.diff(heat) >= 0)
    assert heat[-1] == final['physical_ohmic_heat_j']
    cfg = record['configuration']
    if cfg['resistivity'] == 0:
        assert np.all(heat == 0)
    if cfg.get('dt_max') is not None:
        assert np.max(np.diff(time)) <= cfg['dt_max']*(1+1e-12)
    assert final['max_hall_nonlinear_residual'] <= cfg.get('nonlinear_tolerance', 1e-12)
    result = {}
    for prefix, diagnostics in (('initial', record['initial']), ('final', final)):
        u, psi = archive[prefix+'_conserved'], archive[prefix+'_psi']
        assert u.shape == (mesh.nr, mesh.nz, 6) and psi.shape == (mesh.nr+1, mesh.nz+1)
        assert np.all(u[..., 0] > 0)
        s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
        kinetic = .5*(u[..., 1]**2+u[..., 2]**2+u[..., 4]**2/s)/u[..., 0]
        magnetic = .5*magnetic_product(mesh, psi, u[..., 5], psi, u[..., 5])
        internal = u[..., 3]-kinetic-magnetic
        assert np.all(internal > 0)
        for key, field in (('mass_kg', u[..., 0]), ('total_energy_j', u[..., 3]),
                           ('kinetic_energy_j', kinetic), ('magnetic_energy_j', magnetic),
                           ('internal_energy_j', internal), ('angular_momentum_kg_m2_s', u[..., 4])):
            value, scale = float(np.sum(mesh.volumes*field)), float(np.sum(mesh.volumes*abs(field)))
            error = abs(value-diagnostics[key])/max(scale, np.finfo(float).tiny)
            assert error < 3e-12, (cfg['name'], prefix, key, error)
        np.testing.assert_allclose(psi[0], 0., atol=1e-20, rtol=0)
        np.testing.assert_array_equal(psi[:, -1], psi[:, 0])
        np.testing.assert_allclose(psi[-1], psi[-1, 0], atol=1e-18, rtol=0)
        result[prefix+'_replayed_min_pressure_pa'] = float((gamma-1)*internal.min())
    for key in ('total_energy_j', 'mass_kg'):
        values = np.asarray(history.get(key, [record['initial'][key], final[key]]))
        result['max_'+key+'_relative_drift'] = float(max(abs(values-values[0]))/abs(values[0]))
    angular = np.asarray(history.get('angular_momentum_kg_m2_s', [record['initial']['angular_momentum_kg_m2_s'], final['angular_momentum_kg_m2_s']]))
    result['max_angular_momentum_absolute_drift'] = float(max(abs(angular-angular[0])))
    result['replayed_steps'] = len(time)-1
    return result


def legacy_axis(directory='outputs/hall_work_verification'):
    """Validate historical source against its commit, not current changed files."""
    directory = Path(directory)
    data = json.loads((directory/'measurements.json').read_text())
    for name, sha in data['provenance']['source_files'].items():
        historical = subprocess.check_output(['git', 'show', f'f28714c:{name}'])
        assert hashlib.sha256(historical).hexdigest() == sha
    result = []
    for nr in (8, 16, 32, 64, 128):
        name = f'nonlinear_n{nr}'
        record = next(r for r in data['cases'] if r['configuration']['name'] == name)
        assert record == json.loads((directory/(name+'.json')).read_text())
        assert hashlib.sha256((directory/record['fields']['file']).read_bytes()).hexdigest() == record['fields']['sha256']
        result.append(dict(nr=nr, error=record['metrics']['axis_g_rate_relative_L2']))
    return dict(commit='f28714c', source_sha256=data['provenance']['source_sha256'], first_annulus=result)


def load_all():
    main = collection('outputs/compatible_hall_verification', protocol())
    extra = collection('outputs/compatible_hall_controls', controls())
    return tuple(a | b for a, b in zip(main, extra))
