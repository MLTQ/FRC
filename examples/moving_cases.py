"""Reproducible ideal material-plasma fixtures, independent references and archives."""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import platform
import sys
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma
from plasma_bench.waveforms import Gaussian


@dataclass(frozen=True)
class MovingCase:
    name: str
    kind: str
    radial_cells: int
    end_time: float
    max_timestep: float = .001
    matching_radius: float = 2.
    modes: int = 8
    quadrature: int = 5
    cfl: float = .15

    def build(self):
        n = self.radial_cells
        d = MovingDomain(n, 2*n, matching_radius=self.matching_radius,
                         vacuum_layers=round(n*(self.matching_radius-1)),
                         modes=self.modes, quadrature=self.quadrature)
        common = dict(max_timestep=self.max_timestep, cfl=self.cfl)
        if self.kind == "expansion":
            return MovingPlasma(d, 0., pressure=lambda r,z:1-r*r-z*z, **common)
        if self.kind not in ("frc", "source"):
            raise ValueError("Unknown moving reference")
        bias = np.sqrt(MU0) if self.kind == "frc" else Gaussian(np.sqrt(MU0), .04, .04)
        b = bias if np.isscalar(bias) else float(bias.value(0.))
        return MovingPlasma(d, lambda r,z:-.75*b*(1-r*r-z*z), bias=bias,
                            pressure=lambda r,z:5.625*b*b/MU0*r*r*(1-r*r-z*z), **common)

    def configuration(self):
        return {**asdict(self), "radius_m":1., "angular_cells":2*self.radial_cells,
                "vacuum_layers":round(self.radial_cells*(self.matching_radius-1)),
                "gamma":5/3, "density_kg_m3":1., "initial_velocity_m_s":0.,
                "physical_resistivity_ohm_m":0., "numerical_viscosity":0.,
                "numerical_resistivity":0., "exterior_fluid":False,
                "bias":({"kind":"Gaussian", "amplitude_t":float(np.sqrt(MU0)),
                         "center_s":.04, "sigma_s":.04} if self.kind == "source" else
                        {"kind":"DC", "amplitude_t":float(np.sqrt(MU0)) if self.kind == "frc" else 0.}),
                "initial_a_formula":"0" if self.kind == "expansion" else "-3 B(0) (1-r^2-z^2)/4",
                "initial_pressure_pa":"1-r^2-z^2" if self.kind == "expansion" else
                    "(45/8) B(0)^2/mu0 r^2 (1-r^2-z^2)",
                "interpretation":"Unit-scale numerical fixture; no species/transport closure or formation claim"}


def provenance():
    paths = sorted(Path("plasma_bench").glob("*.py"))+[
        Path("examples/moving_cases.py"), Path("examples/moving_measurements.py"),
        Path("examples/moving_static_measurements.py")]
    hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    combined = hashlib.sha256("\n".join(f"{p}:{h}" for p,h in hashes.items()).encode()).hexdigest()
    return {"source_sha256":combined, "files":hashes, "python":sys.version,
            "platform":platform.platform(), "numpy":np.__version__, "scipy":scipy.__version__}


def snapshot(s):
    e = s.evaluate();d=s.domain
    return {"coordinates":s.state.coordinates.copy(), "velocity":s.state.velocity.copy(),
            "reference_coordinates":d.fluid_reference.copy(), "cells":d.cells[:d.fluid_cell_count].copy(),
            "boundary":d.boundary.copy(), "induced_a":e["magnetic"].induced.copy(),
            "total_material_a":s.flux.coefficients(s.state.coordinates),
            "flux_labels":s.flux.labels.copy(), "node_mass":s.node_mass.copy(),
            "cell_mass":s.cell_mass.copy(), "pressure_invariant":s.pressure_invariant.copy(),
            "density":s.cell_mass/e["volume"], "pressure":e["pressure"],
            "volume":e["volume"], "bias":np.array(e["bias"])}


def errors(s, kind):
    d=s.domain;x0=d.fluid_reference;e=s.evaluate();mass=s.node_mass
    result={"initial_or_current_acceleration_mass_rms_m_s2":float(np.sqrt(
        np.sum(e["force"]**2/mass[:,None])/mass.sum())),
        "nodal_psi_max_drift_wb_per_radian":float(np.max(abs(s.state.coordinates[:,0]**2*
            s.flux.coefficients(s.state.coordinates)-x0[:,0]**2*s.flux.labels)))}
    if kind == "expansion":
        scale=np.sqrt(1+2*s.time**2);speed=2*s.time/scale
        displacement=s.state.coordinates-scale*x0
        dv=s.state.velocity-speed*x0
        result.update(exact_scale=float(scale), position_mass_rms_error_m=float(np.sqrt(
            np.sum(mass[:,None]*displacement**2)/mass.sum())),
            maximum_boundary_position_error_m=float(np.linalg.norm(displacement[d.boundary],axis=-1).max()),
            velocity_relative_mass_L2=float(np.sqrt(np.sum(mass[:,None]*dv**2)/
                np.sum(mass[:,None]*(speed*x0)**2))) if speed else 0.,
            density_relative_volume_L2=float(np.sqrt(np.sum(e["volume"]*(s.cell_mass/e["volume"]*scale**3-1)**2)/e["volume"].sum())),
            volume_relative_error=float(e["volume"].sum()/(4*np.pi/3*scale**3)-1))
    else:
        result.update(maximum_material_displacement_m=float(np.linalg.norm(s.state.coordinates-x0,axis=-1).max()),
            maximum_boundary_displacement_m=float(np.linalg.norm(s.state.coordinates[d.boundary]-x0[d.boundary],axis=-1).max()))
    return result
