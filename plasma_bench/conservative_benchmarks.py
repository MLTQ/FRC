"""Identical initial fixtures for the conservative and split vacuum solvers."""
from dataclasses import asdict
from scipy.constants import mu_0 as MU0
from .conservative_coupled import ConservativeVacuumMHD, ConservativeParameters, ConservativeState


def from_split_initial(solver, **overrides):
    """Convert only initial gas energy to total cell MHD energy; no evolved reset."""
    if solver.time != 0 or solver.steps:
        raise ValueError("Fixture conversion accepts only an unadvanced initial state")
    if solver.parameters.shock_diffusion:
        raise ValueError("Conservative UCT has no separately tuned shock diffusion")
    conserved = solver.state.gas.copy()
    br, bz = solver.mesh.cell_fields(solver.mapping.fluid_flux(solver.total_flux()))
    conserved[..., 3] += (br*br+bz*bz)/(2*MU0)
    options = asdict(solver.parameters)
    options = {key: options[key] for key in ConservativeParameters.__dataclass_fields__ if key in options}
    options.update(overrides)
    params = ConservativeParameters(**options)
    return ConservativeVacuumMHD(solver.mesh, solver.sources, ConservativeState(conserved, solver.state.induced_psi.copy()), params)
