"""Explicit, validated and JSON-serializable experiment configuration."""
from dataclasses import dataclass, asdict, field
import json
from pathlib import Path
import numpy as np
from scipy.constants import atomic_mass
from .geometry import Cylinder, Grid
from .coils import Coil
from .waveforms import Unipolar, to_dict as wave_dict, from_dict as wave_from_dict


@dataclass(frozen=True)
class Material:
    region: Cylinder = field(default_factory=lambda: Cylinder(0.035, 0.10))
    conductivity: float = 2e4
    species: str = "Ar+"
    ion_mass: float = 39.948*atomic_mass
    charge_state: int = 1
    ion_density: float = 1e20
    electron_temperature_ev: float = 5.0
    ion_temperature_ev: float = 1.0
    initial_ionization: float = 1.0

    def __post_init__(self):
        values = (self.conductivity, self.ion_mass, self.ion_density,
                  self.electron_temperature_ev, self.ion_temperature_ev)
        if not all(np.isfinite(v) and v > 0 for v in values):
            raise ValueError("Material parameters must be finite and positive")
        if not isinstance(self.charge_state, int) or self.charge_state < 1:
            raise ValueError("Charge state must be a positive integer")
        if not 0 < self.initial_ionization <= 1:
            raise ValueError("Initial ionization fraction must lie in (0,1]")


def default_coils():
    wave = Unipolar(1000)
    return (Coil(0.061, -0.055, wave, name="formation_left"),
            Coil(0.061, 0.055, wave, name="formation_right"))


@dataclass(frozen=True)
class Experiment:
    grid: Grid = field(default_factory=Grid)
    chamber: Cylinder = field(default_factory=Cylinder)
    material: Material = field(default_factory=Material)
    coils: tuple[Coil, ...] = field(default_factory=default_coils)
    bias: float = 0.005
    duration: float = 100e-6
    dt: float = 0.25e-6
    stage: int = 1
    label: str = "two-coil formation experiment"
    output_stride: int = 4
    persistence_threshold: float = 5e-6

    def __post_init__(self):
        object.__setattr__(self, "coils", tuple(self.coils))
        if not np.isfinite(self.bias) or not all(np.isfinite(v) and v > 0 for v in (self.duration, self.dt)):
            raise ValueError("Bias must be finite; duration and dt must be positive")
        if self.dt > self.duration:
            raise ValueError("dt must not exceed duration")
        if not isinstance(self.output_stride, int) or self.output_stride < 1:
            raise ValueError("output_stride must be a positive integer")
        if not np.isfinite(self.persistence_threshold) or self.persistence_threshold < 0:
            raise ValueError("Persistence threshold must be finite and nonnegative")
        if self.stage not in (0, 1, 2):
            raise ValueError("Only stages 0, 1 and the Stage 2 interface are available")
        p, c, g = self.material.region, self.chamber, self.grid
        if p.radius > c.radius or p.half_length > c.half_length:
            raise ValueError("Plasma must fit inside chamber")
        if c.radius >= g.radius or c.half_length >= g.half_length:
            raise ValueError("Outer computational boundary must be beyond the chamber")
        if p.radius < 3*g.dr or p.half_length < 3*g.dz:
            raise ValueError("Plasma must span at least three grid intervals in each direction")
        for coil in self.coils:
            if coil.radius <= c.radius and abs(coil.z) <= c.half_length:
                raise ValueError("External source filaments must be outside the chamber")

    @property
    def times(self):
        # Avoid a spurious extra step when decimal input makes 400 become 400+1ulp.
        n = int(np.ceil(np.nextafter(self.duration/self.dt, -np.inf)))
        return np.linspace(0, self.duration, n+1)

    def to_dict(self):
        data = asdict(self)
        data["coils"] = [{**{k: v for k, v in asdict(c).items() if k != "waveform"},
                          "waveform": wave_dict(c.waveform)} for c in self.coils]
        return data

    @classmethod
    def from_dict(cls, data):
        d = dict(data)
        if "grid" in d:
            d["grid"] = Grid(**d["grid"])
        if "chamber" in d:
            d["chamber"] = Cylinder(**d["chamber"])
        if "material" in d:
            m = dict(d["material"])
            if "region" in m:
                m["region"] = Cylinder(**m["region"])
            d["material"] = Material(**m)
        if "coils" in d:
            d["coils"] = tuple(Coil(**{**c, "waveform": wave_from_dict(c["waveform"])}) for c in d["coils"])
        return cls(**d)

    def save(self, path):
        Path(path).write_text(json.dumps(self.to_dict(), indent=2)+"\n")

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text()))
