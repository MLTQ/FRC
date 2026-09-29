"""Strict case, source, archive and independent thermodynamic replay audit."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.toroidal_state import ToroidalState
from .warm_reference import WarmReference
from .warm_cases import protocol
from .warm_records import provenance
from .warm_norms import physical_norms
from .compatible_hall_audit import magnetic_product


def history_rows(columns, count):
    """Validate and invert the lossless scalar/nested-column history layout."""
    expanded = {}
    for key, values in columns.items():
        if isinstance(values, dict):
            expanded[key] = history_rows(values, count)
        else:
            assert isinstance(values, list) and len(values) == count
            expanded[key] = values
    return [{key: values[i] for key, values in expanded.items()} for i in range(count)]


def direct_density(mesh, conserved, ion_mass):
    """Explicit integration of incident quarter primal cells, independent loops."""
    n, volume = np.zeros((mesh.nr+1, mesh.nz)), np.zeros((mesh.nr+1, mesh.nz))
    for i in range(mesh.nr+1):
        for j in range(mesh.nz):
            for cell_i in (i-1, i):
                if 0 <= cell_i < mesh.nr:
                    for cell_j in ((j-1) % mesh.nz, j):
                        support = np.pi*(mesh.rf[cell_i+1]**2-mesh.rf[cell_i]**2)*mesh.dz/4
                        volume[i, j] += support
                        n[i, j] += support*conserved[cell_i, cell_j, 0]/ion_mass
    return n/volume, volume


def replay(mesh, reference, scale, arrays, prefix, diagnostic):
    u, psi, entropy = (arrays[prefix+'_'+k] for k in ('conserved', 'psi', 'entropy'))
    n, dual_volume = direct_density(mesh, u, reference.ion_mass)
    pe = (n/scale)**(reference.gamma-1)*entropy
    s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
    kinetic = .5*(u[..., 1]**2+u[..., 2]**2+u[..., 4]**2/s)/u[..., 0]
    magnetic = .5*magnetic_product(mesh, psi, u[..., 5], psi, u[..., 5])
    internal = u[..., 3]-kinetic-magnetic
    pe_cell = np.empty((mesh.nr, mesh.nz))
    for i in range(mesh.nr):
        for j in range(mesh.nz):
            pe_cell[i, j] = sum(pe[a, b] for a in (i, i+1) for b in (j, (j+1) % mesh.nz))/4
    pi = (reference.gamma-1)*internal-pe_cell
    assert np.all(n > 0) and np.all(pe > 0) and np.all(pi > 0)
    np.testing.assert_allclose(dual_volume, arrays['dual_volumes'], rtol=6e-16)
    values = dict(mass_kg=np.sum(mesh.volumes*u[..., 0]), total_energy_j=np.sum(mesh.volumes*u[..., 3]),
                  kinetic_energy_j=np.sum(mesh.volumes*kinetic), magnetic_energy_j=np.sum(mesh.volumes*magnetic),
                  internal_energy_j=np.sum(mesh.volumes*internal), electron_internal_energy_j=np.sum(dual_volume*pe)/(reference.gamma-1),
                  ion_internal_energy_j=np.sum(mesh.volumes*pi)/(reference.gamma-1),
                  electron_entropy_proxy_integral_j=np.sum(dual_volume*entropy),
                  electron_log_entropy_particles=np.sum(dual_volume*n*np.log(entropy/(n*E))))
    for name, value in values.items():
        np.testing.assert_allclose(value, diagnostic[name], rtol=4e-12, atol=1e-28, err_msg=name)
    return pe, pi, n


def audit(directory='outputs/warm_adiabatic_verification', configurations=None):
    directory = Path(directory)
    data = json.loads((directory/'measurements.json').read_text())
    assert data['completed'] and data['provenance'] == provenance()
    assert [r['configuration'] for r in data['cases']] == (protocol() if configurations is None else configurations)
    records, arrays, results = {}, {}, {}
    for record in data['cases']:
        case = record['configuration']; name = case['name']
        assert record['completed'] and record['provenance'] == data['provenance'] and name not in records
        assert record == json.loads((directory/(name+'.json')).read_text())
        path = directory/record['fields']['file']
        assert path.name == name+'.npz' and hashlib.sha256(path.read_bytes()).hexdigest() == record['fields']['sha256']
        with np.load(path, allow_pickle=False) as f:
            archive = {k: f[k] for k in f.files}
        assert all(np.isfinite(v).all() for v in archive.values())
        ref = WarmReference(**case['reference'])
        mesh = MHDMesh(case['radius'], ref.length/2, case['nr'], case['nz'])
        np.testing.assert_array_equal(archive['radial_faces'], mesh.rf)
        np.testing.assert_array_equal(archive['axial_faces'], mesh.zf)
        np.testing.assert_array_equal(archive['primal_volumes'], mesh.volumes)
        initial_exact = ref.cell_means(mesh, 0., 18)
        for key, slot in (('density', 0), ('momentum', 2), ('energy', 3)):
            np.testing.assert_allclose(archive['initial_conserved'][..., slot], initial_exact[key], rtol=1e-13, atol=1e-20)
        replay(mesh, ref, case['density_scale'], archive, 'initial', record['initial'])
        pe, pi, n = replay(mesh, ref, case['density_scale'], archive, 'final', record['final'])
        np.testing.assert_allclose(pe, archive['final_electron_pressure'], rtol=4e-15)
        np.testing.assert_allclose(pi, archive['final_ion_pressure'], rtol=4e-13)
        history = history_rows(record['history_columns'], record['accepted_steps']+1)
        time = np.array([h['time_s'] for h in history])
        assert history[0] == record['initial'] and history[-1] == record['final']
        assert len(history) == record['accepted_steps']+1 and np.all(np.diff(time) > 0) and time[0] == 0
        assert abs(time[-1]-record['end_time_s']) < 2e-13*record['end_time_s']
        assert time[-1] == record['actual_time_s']
        for i, h in enumerate(history):
            assert h['steps'] == i and h['physical_ohmic_heat_j'] == 0
            assert all(h[k] > 0 for k in ('min_density_kg_m3', 'min_pressure_pa', 'min_electron_pressure_pa', 'min_ion_pressure_pa'))
            assert abs(h['electron_internal_energy_j']+h['ion_internal_energy_j']-h['internal_energy_j']) < 3e-15*h['total_energy_j']
        if case['fixed_steps'] is not None:
            assert len(time)-1 == case['fixed_steps']
            np.testing.assert_allclose(np.diff(time), record['end_time_s']/case['fixed_steps'], rtol=5e-13)
        state = ToroidalState(archive['final_conserved'], archive['final_psi'], 'compatible')
        norms = physical_norms(mesh, state, pe, n, ref, time[-1], reconstruction_order=case['order'])
        fine = physical_norms(mesh, state, pe, n, ref, time[-1], order=12, reconstruction_order=case['order'])
        result = dict(physical_norms=norms, quadrature_relative_norm_change=max(abs(norms[k]['perturbation_relative_L2']-fine[k]['perturbation_relative_L2']) for k in norms))
        for key in ('mass_kg', 'total_energy_j', 'electron_entropy_proxy_integral_j'):
            v = np.array([h[key] for h in history])
            result['max_'+key+'_relative_drift'] = float(np.max(abs(v-v[0]))/abs(v[0]))
        result['max_log_entropy_change_per_particle'] = float(max(abs(h['electron_log_entropy_particles_change']) for h in history)/(history[0]['mass_kg']/ref.ion_mass))
        records[name], arrays[name], results[name] = record, archive, result
    return data, records, arrays, results
