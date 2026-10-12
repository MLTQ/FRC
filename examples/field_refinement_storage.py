"""Source-derived retained-storage estimates for separate field refinements."""
import gc
import json
from pathlib import Path
import resource
import time
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import block_diag, csr_matrix, eye, kron
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.mhd_mesh import MHDMesh
from .regular_plasma_reference import reference
from .implicit_plasma_records import provenance, freeze


def sparse_record(matrix):
    return dict(shape=list(matrix.shape), nnz=matrix.nnz,
        bytes=matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)


def measure(nr, nz=32):
    ref = reference(dict(kind='contact'))
    geo = StrongWallGeometry(CompatibleGeometry(MHDMesh(ref.radius, .5*ref.length, nr, nz), quadrature_order=12))
    nodes = (nr+1)*nz
    active = np.r_[np.arange(nz, nr*nz), np.arange(nodes, nodes+nr*nz)]
    radial = null_space(geo.wall_trace[None, :])
    basis = block_diag((eye((nr-1)*nz), kron(csr_matrix(radial), eye(nz))), format='csr')
    full = csr_matrix((np.ones(len(active)), (active, np.arange(len(active)))), shape=(geo.hessian.shape[0], len(active)))@basis
    metric = (full.T@geo.hessian@full).tocsc()
    lift = full.multiply((1/np.sqrt(metric.diagonal()))[None, :]).tocsr()
    dual = (geo.hessian@lift).tocsr()
    matrices = dict(raw_curl=geo.original.faraday, projected_curl=geo.faraday,
        projection=geo.projection, hessian=geo.hessian, lift=lift, dual_lift=dual)
    matrices.update({f'electric_mass_{i}': m for i, m in enumerate(geo.mass)})
    retained = {k: sparse_record(v) for k, v in matrices.items()}
    factors = {f'mass_LU_{i}': dict(L=sparse_record(solve.__self__.L), U=sparse_record(solve.__self__.U))
        for i, solve in enumerate(geo.solves)}
    return dict(radial_cells=nr, probe_axial_cells=nz, matrices=retained,
        probe_mass_factorization=factors, probe_scalar_quadrature_bytes=geo.weight.nbytes,
        process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    identity = provenance(['examples.field_refinement_storage'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=False, source_sha256=identity['source_sha256'], probes=[], estimates=[])
    def save(): (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    save(); start = time.monotonic()
    for nr in (256, 512):
        print('Source geometry storage probe radial', nr, 'axial32', flush=True)
        result['probes'].append(measure(nr)); gc.collect(); save()
    for nr, nz in [(256, 256), (256, 512), (512, 256)]:
        probe = next(p for p in result['probes'] if p['radial_cells'] == nr)
        matrix_bytes = {k: int(v['bytes']*nz/probe['probe_axial_cells']) for k, v in probe['matrices'].items()}
        # UniformWallMetric stores all nz complex128 Cholesky blocks.
        modal_bytes = nz*(2*nr-2)**2*16
        scalar_bytes = nr*nz*12**2*8
        result['estimates'].append(dict(radial_cells=nr, axial_cells=nz,
            projected_matrix_bytes_estimated=matrix_bytes,
            stored_uniform_metric_bytes=modal_bytes,
            scalar_quadrature_bytes=scalar_bytes, vector_quadrature_bytes=3*scalar_bytes,
            vector_gradient_quadrature_bytes=9*scalar_bytes,
            estimated_sparse_plus_metric_GiB=(sum(matrix_bytes.values())+modal_bytes)/2**30))
    if provenance(['examples.field_refinement_storage']) != identity: raise ValueError('Source changed')
    result.update(terminal=True, elapsed_seconds=time.monotonic()-start,
        scope='Storage accounting, not a physics solve or peak-memory guarantee. Sparse extrapolation assumes fixed axial stencil; large-grid LU fill, fluid/current/species temporaries and allocator workspaces excluded.',
        rss_units='macOS ru_maxrss bytes; cumulative process maximum, not separate per-case peak')
    save(); return result
