"""Declared wall-conforming initialization of independent continuum fixtures."""
import numpy as np
from plasma_bench.regular_field_map import RegularFieldMap
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from plasma_bench.strong_wall_geometry import StrongWallGeometry
from plasma_bench.strong_wall_plasma import StrongWallFluid, StrongWallPlasma
from .regular_plasma_reference import initialize as original_initialize


def initialize(spec):
    ref, original_fluid, u, original, field = original_initialize(spec)
    fluid = StrongWallFluid(ref.radius, ref.length, original_fluid.elements,
        original_fluid.nz, original_fluid.mass, degree=original_fluid.degree,
        gammas=tuple(original_fluid.thermal.gammas))
    # Analytic J1(qR)=0 has floating-point residue. This inactive slot belongs
    # to the initial approximation, not to an accepted-state modification.
    u[-1, -1, :, 1] = 0.
    geometry = StrongWallGeometry(original.geometry)
    system = StrongWallPlasma(RegularFieldMap(fluid, geometry), ref.ion_mass,
        **original.solver_options)
    system.original_system = original
    g = geometry.radial_projection@field.g
    field = HarmonicMagnetic(field.mesh, field.psi_increment, g, field.axial_bias)
    return ref, fluid, u, system, field
