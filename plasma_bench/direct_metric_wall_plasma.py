"""Identical split-inverse wall solve without an unused intermediate LU."""
from copy import copy
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import block_diag, csr_matrix, eye, kron
from .split_inverse_wall_induction import SplitInverseWallInduction
from .split_inverse_wall_plasma import SplitInverseWallPlasma
from .uniform_wall_metric import UniformWallMetric
from .current_transport import CurrentTransport
from .regular_species import RegularSpecies


class DirectMetricWallInduction(SplitInverseWallInduction):
    def __init__(self, model, *, tolerance=2e-12, max_iterations=1000):
        if not np.isfinite(tolerance) or not 0 < tolerance < 1 or type(max_iterations) is not int or max_iterations < 1:
            raise ValueError('Positive finite tolerance below one and iteration budget required')
        self.model = model; self.tolerance = tolerance; self.max_iterations = max_iterations
        geo = model.geometry; mesh = model.mesh; nodes = (mesh.nr+1)*mesh.nz
        self.active = np.r_[np.arange(mesh.nz, mesh.nr*mesh.nz), np.arange(nodes, nodes+mesh.nr*mesh.nz)]
        radial = null_space(geo.wall_trace[None, :])
        basis = block_diag((eye((mesh.nr-1)*mesh.nz), kron(csr_matrix(radial), eye(mesh.nz))), format='csr')
        full = csr_matrix((np.ones(len(self.active)), (self.active, np.arange(len(self.active)))),
                          shape=(geo.hessian.shape[0], len(self.active)))@basis
        metric = (full.T@geo.hessian@full).tocsc()
        scale = 1/np.sqrt(metric.diagonal())
        self.lift = full.multiply(scale[None, :]).tocsr()
        self.dual_lift = (geo.hessian@self.lift).tocsr()
        # UniformWallInduction immediately replaces the base LU with this metric.
        # Omit only the unused factorization; all coordinates and solve code remain.
        self.factor = UniformWallMetric(model, self.lift)


class DirectMetricWallPlasma(SplitInverseWallPlasma):
    def operators(self, model):
        if self._operators is None or not self.cache_geometry:
            operators = (CurrentTransport(model), RegularSpecies(model), DirectMetricWallInduction(model, **self.solver_options))
            if not self.cache_geometry:
                return operators
            self._operators = operators
        operators = tuple(copy(operator) for operator in self._operators)
        for operator in operators:
            operator.model = model
        return operators


def direct_metric(system):
    result = DirectMetricWallPlasma(system.mapping, system.ion_mass,
        cache_geometry=system.cache_geometry, **system.solver_options)
    result.original_system = system.original_system
    return result
