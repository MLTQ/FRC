"""Rebuild and compare full curl actions at every archived rejected iterate."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import block_diag, csr_matrix, eye, kron
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.mhd_mesh import MHDMesh
from .regular_plasma_reference import reference
from .curl_evaluation_actions import evaluate
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, original='outputs/split_residual_probe'):
    root = Path(original); destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    record = json.loads((root/'summary.json').read_text())
    source = json.loads((root/'provenance.json').read_text())
    if not record['terminal'] or not record['observation_consistent']:
        raise ValueError('Complete original matched rejection required')
    if record['source_sha256'] != source['source_sha256'] or any(sha(p) != h for p, h in source['source_files'].items()):
        raise ValueError('Unchanged original source closure required')
    path = root/record['arrays_file']
    if sha(path) != record['arrays_sha256']: raise ValueError('Unchanged full observations required')
    saved = dict(np.load(path)); identity = provenance(['examples.curl_evaluation_probe'])
    destination.mkdir(parents=True); freeze(destination, identity); start = time.monotonic()
    result = dict(terminal=False, source_sha256=identity['source_sha256'], input_directory=str(root),
        input_summary_sha256=sha(root/'summary.json'), input_arrays_sha256=sha(path), observations=[])
    def save(): (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    save(); print('Rebuilding original full geometry and dual lift', flush=True)
    spec = record['spec']; ref = reference(dict(spec, kind='contact'))
    mesh = MHDMesh(ref.radius, .5*ref.length, spec['field_nr'], spec['field_nz'])
    geo = StrongWallGeometry(CompatibleGeometry(mesh, degree=5, quadrature_order=spec.get('quadrature', 12)))
    nodes = (mesh.nr+1)*mesh.nz
    active = np.r_[np.arange(mesh.nz, mesh.nr*mesh.nz), np.arange(nodes, nodes+mesh.nr*mesh.nz)]
    radial = null_space(geo.wall_trace[None, :])
    basis = block_diag((eye((mesh.nr-1)*mesh.nz), kron(csr_matrix(radial), eye(mesh.nz))), format='csr')
    full = csr_matrix((np.ones(len(active)), (active, np.arange(len(active)))), shape=(geo.hessian.shape[0], len(active)))@basis
    metric = (full.T@geo.hessian@full).tocsc()
    lift = full.multiply((1/np.sqrt(metric.diagonal()))[None, :]).tocsr()
    dual = (geo.hessian@lift).tocsr()
    del full, basis, metric
    arrays = {}; norm = np.linalg.norm
    for index in range(len(record['trace'])):
        print('Comparing complete actual curl actions, refinement', index, flush=True)
        prefix = f'refinement_{index}_'; data = {k[len(prefix):]: v for k, v in saved.items() if k.startswith(prefix)}
        actions = {name: evaluate(geo, data[key], exact=index == len(record['trace'])-1)
            for name, key in [('e0', 'forcing_electric'), ('ei', 'inertia_electric'), ('e', 'electric'), ('low', 'electric_addition_low')]}
        rhs = dual.T@actions['e0']['assembled']; scale = max(norm(rhs), np.finfo(float).tiny)
        physical = dual.T@(actions['e']['assembled']-data['magnetic_rate'])
        distribution = dual.T@(actions['e']['assembled']-actions['e0']['assembled']-actions['ei']['assembled'])
        identities = dict(rhs=np.array_equal(rhs, data['rhs']),
            physical_dual=np.array_equal(physical, data['physical_dual_residual']),
            electric_distribution=np.array_equal(distribution, data['electric_curl_distribution']),
            magnetic_rate=np.array_equal(lift@data['solution'], data['magnetic_rate']))
        if not all(identities.values()): raise ValueError('Rebuilt actual archived actions differ')
        metrics = {}
        for method in actions['e']:
            if method in ('raw', 'paired_raw', 'exact_raw'): continue
            f0, fi, fe, fl = [actions[name][method] for name in ('e0', 'ei', 'e', 'low')]
            gap = dual.T@(fe-f0-fi); difference = fe-data['magnetic_rate']
            fieldnorm = np.sqrt(data['magnetic_rate']@(geo.hessian@data['magnetic_rate']))
            metrics[method] = dict(distribution_over_rhs=float(norm(gap)/scale),
                distribution_plus_low_over_rhs=float(norm(gap+dual.T@fl)/scale),
                full_action_change_over_rhs=float(norm(dual.T@(fe-actions['e']['assembled']))/scale),
                forcing_action_change_over_rhs=float(norm(dual.T@(f0-actions['e0']['assembled']))/scale),
                dual_faraday_original_rhs=float(norm(dual.T@difference)/scale),
                physical_faraday=float(np.sqrt(difference@(geo.hessian@difference))/max(fieldnorm, np.finfo(float).tiny)))
            arrays[prefix+method+'_distribution'] = gap
            arrays[prefix+method+'_faraday_difference'] = difference
        exact_metrics = {}
        if 'exact_raw' in actions['e']:
            for name, values in actions.items():
                exact_metrics[name] = {key: float(norm(dual.T@geo.project_curl(values[key]-values['exact_raw']))/scale)
                    for key in ('raw', 'paired_raw')}
        for name, values in actions.items():
            arrays.update({prefix+name+'_'+key: value for key, value in values.items()})
        result['observations'].append(dict(refinement=index, original_bitwise_identities=identities,
            metrics=metrics, raw_error_against_exact_over_rhs=exact_metrics))
        save()
    path = destination/'actions.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.curl_evaluation_probe']) != identity: raise ValueError('Source changed')
    result.update(terminal=True, arrays_file=path.name, arrays_sha256=sha(path),
        elapsed_seconds=time.monotonic()-start, full_plasma_rate_accepted=False,
        exact_scope='Final iterate full raw sparse products/sums only; projected actions remain binary64',
        scope='Actual rejected iterates; alternate residuals do not accept or repair a full plasma derivative')
    save(); print('Complete curl evaluation observation terminal', flush=True)
    return result
