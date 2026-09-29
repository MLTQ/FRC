"""Frozen provenance and complete moving-resistive experiment specifications."""
from dataclasses import dataclass,asdict
from pathlib import Path
import hashlib
import platform
import sys
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.resistive_moving import ResistiveMovingPlasma
from plasma_bench.material_resistivity import MaterialResistivity
from plasma_bench.waveforms import Gaussian
from .moving_cases import snapshot as ideal_snapshot


def provenance():
    names=("resistive_records","resistive_decay_measurements","resistive_dynamic_measurements",
           "spherical_decay_reference","moving_cases")
    paths=sorted(Path("plasma_bench").glob("*.py"))+[Path(f"examples/{n}.py") for n in names]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    digest=hashlib.sha256("\n".join(f"{p}:{v}" for p,v in hashes.items()).encode()).hexdigest()
    return {"source_sha256":digest,"files":hashes,"python":sys.version,"platform":platform.platform(),
            "numpy":np.__version__,"scipy":scipy.__version__}


@dataclass(frozen=True)
class ResistiveCase:
    name: str
    kind: str
    radial_cells: int
    diffusivity: float = .1
    end_time: float = .04
    max_timestep: float = .0005
    diffusion_tolerance: float = 1e-4
    matching_radius: float = 2.
    modes: int = 8
    quadrature: int = 5
    cfl: float = .15

    def build(self):
        if self.kind not in ("frc","source"):
            raise ValueError("Unknown moving resistive fixture")
        n=self.radial_cells
        d=MovingDomain(n,2*n,matching_radius=self.matching_radius,
                       vacuum_layers=round(n*(self.matching_radius-1)),modes=self.modes,quadrature=self.quadrature)
        bias=np.sqrt(MU0) if self.kind=="frc" else Gaussian(np.sqrt(MU0),.04,.04)
        b=float(bias if np.isscalar(bias) else bias.value(0.))
        return ResistiveMovingPlasma(d,lambda r,z:-.75*b*(1-r*r-z*z),
            pressure=lambda r,z:5.625*b*b/MU0*r*r*(1-r*r-z*z),bias=bias,
            resistivity=self.diffusivity*MU0,diffusion_tolerance=self.diffusion_tolerance,
            cfl=self.cfl,max_timestep=self.max_timestep)

    def configuration(self):
        return {**asdict(self),"angular_cells":2*self.radial_cells,"radius_m":1.,
                "vacuum_layers":round(self.radial_cells*(self.matching_radius-1)),"density_kg_m3":1.,
                "gamma":5/3,"resistivity_ohm_m":self.diffusivity*MU0,
                "initial_velocity_m_s":0.,"initial_a":"-3 B(0) (1-r^2-z^2)/4",
                "initial_pressure_pa":"(45/8) B(0)^2/mu0 r^2 (1-r^2-z^2)",
                "bias":({"kind":"Gaussian","amplitude_t":float(np.sqrt(MU0)),"center_s":.04,"width_s":.04}
                        if self.kind=="source" else {"kind":"DC","amplitude_t":float(np.sqrt(MU0))}),
                "exterior_fluid":False,"numerical_viscosity":0.,"numerical_resistivity":0.,
                "interpretation":"Unit-scale scalar-eta moving verification; no species/transport/Hall closure or formation claim"}


def snapshot(s):
    result=ideal_snapshot(s)
    operator=MaterialResistivity(s.domain,s.state.coordinates,s.resistivity)
    q=operator.current_over_radius(result["total_material_a"],float(s.bias.value(s.time)))
    result.update(current_over_radius=q,current_quadrature=operator.current_quadrature(q),
                  quadrature_coordinates=operator.elements[0].points[:s.domain.fluid_cell_count],
                  cell_ohmic_heat=s.cell_ohmic_heat.copy())
    return result
