"""Diagonal equilibration of positive physical mass systems, without regularization."""
import numpy as np
from scipy.linalg import eigh, solve_banded


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


def largest_mass_eigenvalue(stiffness, mass):
    """Largest eigenvalue of K v=lambda M v in balanced physical coordinates."""
    stiffness, mass = np.asarray(stiffness, float), np.asarray(mass, float)
    if (mass.ndim != 2 or mass.shape[0] != mass.shape[1] or not len(mass)
            or stiffness.shape != mass.shape or not np.isfinite(mass).all()
            or not np.isfinite(stiffness).all() or np.any(np.diag(mass) <= 0)):
        raise ValueError('Finite square physical matrices with positive mass diagonal required')
    scale = np.sqrt(np.diag(mass))
    balanced_mass = (mass/scale[:, None])/scale[None, :]
    balanced_stiffness = (stiffness/scale[:, None])/scale[None, :]
    if not np.isfinite(balanced_mass).all() or not np.isfinite(balanced_stiffness).all():
        raise ValueError('Balanced physical eigenproblem exceeds floating-point range')
    magnitude = float(abs(balanced_stiffness).max())
    if magnitude == 0:
        return 0.
    value = eigh(balanced_stiffness/magnitude, balanced_mass,
        subset_by_index=[len(mass)-1, len(mass)-1], eigvals_only=True,
        driver='gvx', check_finite=False)[0]*magnitude
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Nonzero physical stiffness requires a representable positive largest eigenvalue')
    return float(value)
