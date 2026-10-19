"""A posteriori local conservative/CT fallback with a checked baseline."""
from dataclasses import dataclass
import numpy as np
from .mhd_state import NonphysicalState
from .magnetic_quadrature import selected_energy


@dataclass
class LimiterReport:
    troubled: np.ndarray
    affected: np.ndarray
    iterations: int = 0
    minimum_weight: float = 1.


def invalid_cells(mesh, conserved, psi, energy_quadrature="cell", exempt=None):
    """Strict positive density and gas internal energy, with no absolute floor.

    `exempt` cells (dual-energy cold cells) need positive density only; their
    pressure comes from the separately advected entropy.
    """
    rho = conserved[..., 0]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        vr, vz = conserved[..., 1]/rho, conserved[..., 2]/rho
        internal = conserved[..., 3]-.5*rho*(vr*vr+vz*vz)-selected_energy(mesh, psi, energy_quadrature)
    thermal_bad = ~np.isfinite(internal) | (internal <= 0)
    if exempt is not None:
        thermal_bad &= ~exempt
    return ~np.all(np.isfinite(conserved), axis=-1) | thermal_bad | (rho <= 0)


def admissible_update(mesh, conserved, psi, dt, high, low_factory, *, mode="local", energy_quadrature="cell",
                      exempt=None):
    """Return an admissible RHS and report; call baseline factory only if needed.

    Off checks but never limits. Global retains the earlier convex blend as an
    explicit comparison. Local replacement terminates conditionally on an
    admissible first-order proposal, since each pass adds a new troubled cell.
    """
    if mode not in ("local", "global", "off"):
        raise ValueError("Flux limiter must be local, global or off")
    high_rhs = high.rhs(mesh)
    def bad(rhs):
        return invalid_cells(mesh, conserved+dt*rhs[0], psi+dt*rhs[1], energy_quadrature, exempt)
    troubled = bad(high_rhs)
    if not np.any(troubled):
        return high_rhs, LimiterReport(troubled, troubled.copy())
    if mode == "off":
        raise NonphysicalState("Conservative Euler proposal is nonphysical; no flux fallback enabled")
    low = low_factory()
    low_rhs = low.rhs(mesh)
    if np.any(bad(low_rhs)):
        raise NonphysicalState("First-order conservative Euler proposal is nonphysical; no state repair applied")
    if mode == "global":
        blend = lambda theta: tuple(a+theta*(b-a) for a, b in zip(low_rhs, high_rhs))
        lower, upper = 0., 1.
        for _ in range(40):
            middle = .5*(lower+upper)
            if np.any(bad(blend(middle))):
                upper = middle
            else:
                lower = middle
        theta = .9*lower
        return blend(theta), LimiterReport(troubled, np.ones_like(troubled), 40, theta)
    iterations = 0
    while True:
        mixed, affected = high.local_fallback(low, troubled)
        rhs = mixed.rhs(mesh)
        invalid = bad(rhs)
        iterations += 1
        if not np.any(invalid):
            return rhs, LimiterReport(troubled, affected, iterations, 0.)
        # Flagged cells use the independently checked low update exactly.
        if np.any(invalid & troubled):
            raise NonphysicalState("Local fallback disagrees with its admissible baseline")
        troubled |= invalid
