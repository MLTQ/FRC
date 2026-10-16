"""Diagonal equilibration of positive physical mass systems, without regularization."""
import numpy as np
from scipy.linalg import solve_banded


def solve_mass_banded(matrix, load):
    """Solve M v = b using D=sqrt(diag(M)), A=D^-1 M D^-1.

    The input is a symmetric tridiagonal mass matrix in SciPy band storage.
    No mass, density, eigenvalue or right-hand-side cutoff is introduced.
    """
    matrix, load = np.asarray(matrix, float), np.asarray(load, float)
    if (matrix.ndim != 2 or matrix.shape[0] != 3 or
            load.ndim not in (1, 2) or load.shape[0] != matrix.shape[1] or
            not np.isfinite(matrix).all() or not np.isfinite(load).all() or
            np.any(matrix[1] <= 0)):
        raise ValueError('Finite physical mass/load and positive mass diagonal required')
    scale = np.sqrt(matrix[1])
    balanced = matrix.copy()
    balanced[1] = 1.
    balanced[0, 1:] = matrix[0, 1:]/scale[:-1]/scale[1:]
    balanced[2, :-1] = matrix[2, :-1]/scale[1:]/scale[:-1]
    divisor = scale if load.ndim == 1 else scale[:, None]
    result = solve_banded((1, 1), balanced, load/divisor, check_finite=True)/divisor
    if not np.isfinite(result).all():
        raise ValueError('Physical mass solution exceeds floating-point range')
    return result
