"""Alternative actions and exact raw-curl reference on actual electric arrays."""
from decimal import Decimal, localcontext
import numpy as np


def exact_sparse_action(matrix, value):
    """Correctly round exact products/sums of finite binary64 CSR operands."""
    matrix = matrix.tocsr(); value = np.asarray(value, float)
    if not np.isfinite(value).all() or not np.isfinite(matrix.data).all():
        raise ValueError('Finite binary64 operands required')
    result = np.empty(matrix.shape[0])
    with localcontext() as context:
        context.prec = 2200
        coefficients = [Decimal.from_float(float(x)) for x in matrix.data]
        vectors = [Decimal.from_float(float(x)) for x in value]
        for row, (start, stop) in enumerate(zip(matrix.indptr[:-1], matrix.indptr[1:])):
            result[row] = float(sum((coefficients[k]*vectors[matrix.indices[k]]
                for k in range(start, stop)), Decimal(0)))
    return result


def evaluate(geometry, value, *, exact=False):
    """No changed state: return each separately evaluated full curl vector."""
    raw = geometry.original.faraday@value
    results = dict(assembled=geometry.faraday@value,
        sparse_factored=geometry.projection@raw,
        dense_factored=geometry.project_curl(raw), raw=raw)
    mesh = geometry.mesh; nodes = (mesh.nr+1)*mesh.nz
    er, et, ez = geometry.split_electric(value)
    er = er.reshape(mesh.nr, mesh.nz); ez = ez.reshape(mesh.nr+1, mesh.nz)
    paired = np.zeros_like(raw)
    paired[mesh.nz:mesh.nr*mesh.nz] = (-mesh.rf[1:-1, None]*et.reshape(mesh.nr-1, mesh.nz)).ravel()
    paired[nodes:] = ((er-np.roll(er, -1, axis=1))*(1/mesh.dz)+
        np.diff(ez, axis=0)*(1/mesh.radial_measure[:, None])).ravel()
    results.update(paired_raw=paired, paired_factored=geometry.project_curl(paired))
    if exact:
        results['exact_raw'] = exact_sparse_action(geometry.original.faraday, value)
        results['exact_raw_factored'] = geometry.project_curl(results['exact_raw'])
    return results
