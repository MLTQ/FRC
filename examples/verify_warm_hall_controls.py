"""Audit the supplemental successful controls and preserve all rejected runs."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd_mesh import MHDMesh
from .warm_hall_audit import audit, isolated_endpoints
from .warm_hall_controls import protocol, source_hash
from .warm_hall_records import provenance
from .warm_hall_comparisons import difference
from .verification_records import write_record


def rejected_audit(directory):
    results, records = {}, []
    original = next(c for c in protocol() if c['name'] == 'pressure_fine_mesh_128')
    driver = hashlib.sha256(Path('examples/warm_pressure_failure.py').read_bytes()).hexdigest()
    for count in (320, 640, 1280):
        name = f'pressure_mesh128_time{count}'; record = json.loads((directory/(name+'.json')).read_text())
        case = record['configuration']; ref = case['reference']
        assert case == original | dict(name=name, steps=count) and record['original_control_configuration'] == original
        assert record['experiment_complete'] and not record['completed'] and record['rejected_trial'] is not None
        assert record['provenance'] == provenance() and record['driver_sha256'] == driver and record['control_source_sha256'] == source_hash()
        path = directory/record['fields']['file']; assert path.name == name+'.npz'
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['fields']['sha256']
        with np.load(path, allow_pickle=False) as opened:
            a = {k: opened[k] for k in opened.files}
        assert all(np.isfinite(v).all() for v in a.values())
        mesh = MHDMesh(ref['radius'], ref['length']/2, case['nr'], case['nz'])
        np.testing.assert_array_equal(a['radial_faces'], mesh.rf)
        np.testing.assert_array_equal(a['axial_faces'], mesh.zf)
        np.testing.assert_array_equal(a['primal_volumes'], mesh.volumes)
        isolated_endpoints(record, a, mesh)
        h = record['history_columns']; time = np.array(h['time_s']); dt = record['end_time_s']/count
        assert all(len(v) == record['accepted_steps']+1 and np.isfinite(v).all() for v in h.values())
        for key, v in h.items():
            assert v[0] == record['initial'][key] and v[-1] == record['final'][key]
        np.testing.assert_allclose(time, dt*np.arange(len(time)), rtol=2e-15, atol=0.)
        assert time[-1] == record['actual_time_s'] < record['end_time_s']
        assert min(h['min_ion_pressure_pa']) > 0 and min(h['min_electron_pressure_pa']) > 0
        assert record['rejected_trial']['time_s'] == dt*len(time)
        pe = (a['number_density']/1e20)**(ref['gamma']-1)*a['rejected_entropy']
        pcell = .25*(pe[:-1]+pe[1:]+np.roll(pe[:-1]+pe[1:], -1, axis=1))
        s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
        pi = (ref['gamma']-1)*(a['rejected_total_energy']-.5*s*a['rejected_g']**2/MU0)-pcell
        np.testing.assert_allclose(pi, a['rejected_ion_pressure'], rtol=2e-12, atol=1e-12)
        assert pi.min() < 0 < pe.min() and record['rejected_trial']['nonlinear_residual'] <= case['tolerance']
        np.testing.assert_allclose(pi.min(), record['rejected_trial']['min_ion_pressure_pa'], rtol=2e-12, atol=1e-12)
        energy = np.array(h['total_energy_j'])
        results[name] = dict(accepted_steps=record['accepted_steps'], last_accepted_time_s=time[-1],
                             rejected_trial=record['rejected_trial'], max_total_energy_relative_drift=float(np.max(abs(energy-energy[0]))/energy[0]))
        records.append(record)
    symbol = json.loads((directory/'principal_symbol.json').read_text())
    assert symbol['reference'] == original['reference']
    for path, sha in symbol['sources'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha
    eig = np.linalg.eigvals(1j*np.array(symbol['real_principal_matrix']))
    np.testing.assert_allclose(max(eig.real), symbol['growth_rate_per_wavenumber_m_s'], rtol=1e-14)
    return results, records, symbol


def run():
    _, base, arrays, meshes, _ = audit()
    cases = protocol(); successful = [c for c in cases if c['name'] != 'pressure_fine_mesh_128']
    _, extra, ea, em, audits = audit('outputs/warm_hall_controls', successful, assembled=False)
    assert all(r['control_source_sha256'] == source_hash() for r in extra.values())
    directory = Path('outputs/warm_pressure_failure')
    rejected, failures, symbol = rejected_audit(directory)
    assert {c['name'] for c in cases} == set(extra) | {failures[0]['original_control_configuration']['name']}
    base.update(extra); arrays.update(ea); meshes.update(em)
    times = [difference(base, arrays, meshes, f'pressure_fine_time_{n}', 'pressure_fine_time_640') for n in (160, 320)]
    half = difference(base, arrays, meshes, 'slow_wave_half_cfl', 'wave_b0_16')
    ref = base['wave_b0_16']['configuration']['reference']; pe0 = ref['electron_fraction']*ref['pressure']
    a, b = arrays['slow_wave_amplitude_01']['final_electron_pressure'], arrays['wave_b0_16']['final_electron_pressure']
    volume = arrays['wave_b0_16']['dual_volumes']
    amplitude = float(np.sqrt(np.sum(volume*((a-pe0)/.1-(b-pe0))**2)/np.sum(volume*(b-pe0)**2)))
    from .warm_failure_figure import figure
    path = figure(directory, failures, base)
    sources = ['examples/'+n+'.py' for n in ('verify_warm_hall_controls', 'warm_failure_figure', 'warm_pressure_failure',
               'pressure_split_symbol', 'warm_hall_controls', 'warm_hall_audit', 'warm_hall_field_audit',
               'warm_hall_comparisons', 'warm_audit', 'compatible_hall_audit')]
    report = dict(completed=True, provenance=provenance(), successful_control_audits=audits, rejected_replays=rejected,
                  pressure_fine_time_comparisons=times, slow_half_cfl=half, slow_amplitude_nodal_relative_difference=amplitude,
                  pressure_split_principal_symbol=symbol, figure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  sources={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sources},
                  scope='Five completed controls and three rejected replays, including the original failed sixth control; no formation validation')
    write_record(directory/'report.json', report)
    print('Audited five controls, three rejections, principal symbol. Slow amplitude difference', amplitude, flush=True)


if __name__ == '__main__':
    run()
