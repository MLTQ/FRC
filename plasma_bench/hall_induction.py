"""Experimental nonlinear Hall induction with an exact assembled work identity."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E


@dataclass(frozen=True)
class HallInductionResult:
    psi_rate: np.ndarray
    g_rate: np.ndarray
    electric_coefficients: np.ndarray
    current_samples: np.ndarray
    electric_samples: np.ndarray
    magnetic_work_w: float
    absolute_work_scale_w: float


def hall_induction(geometry, psi, g, number_density):
    """Project J cross B/(ne e); no heat/velocity update and no carrier floor.

    Density may be positive scalar, cell values, or values at geometry's
    physical quadrature. A cell array is explicitly piecewise constant.
    """
    mesh = geometry.mesh
    psi, g, density = np.asarray(psi), np.asarray(g), np.asarray(number_density)
    if psi.shape != (mesh.nr+1, mesh.nz+1) or g.shape != (mesh.nr, mesh.nz) or not np.isfinite(psi).all() or not np.isfinite(g).all():
        raise ValueError('Finite corner psi and cell g required')
    scale = max(float(np.max(abs(psi))), np.finfo(float).tiny)
    if np.max(abs(psi[0])) > 1e-13*scale or np.ptp(psi[-1]) > 1e-13*scale or np.max(abs(psi[:, 0]-psi[:, -1])) > 1e-13*scale:
        raise ValueError('Regular axis, constant wall flux and periodic psi required')
    if density.shape == g.shape:
        density = density[..., None, None]
    density = np.broadcast_to(density, geometry.r.shape)
    if not np.isfinite(density).all() or np.any(density <= 0):
        raise ValueError('Strictly positive carrier density inside this filled fixture required; no floor applied')
    current = geometry.current(psi, g)
    magnetic = geometry.magnetic_samples(psi, g)
    electric = np.cross(current, magnetic)/(density[..., None]*E)
    coefficients = geometry.project_electric(electric)
    rate = geometry.faraday @ coefficients
    psi_rate, g_rate = geometry.unpack(rate)
    work = float(np.dot(geometry.magnetic_gradient(psi, g), rate))
    scale = float(np.sum(geometry.weight*np.linalg.norm(current, axis=-1)*np.linalg.norm(electric, axis=-1)))
    return HallInductionResult(psi_rate, g_rate, coefficients, current, electric, work, scale)
