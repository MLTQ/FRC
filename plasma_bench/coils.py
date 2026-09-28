"""Axisymmetric filamentary source coils, independent of the plasma solver."""
from dataclasses import dataclass
import numpy as np
from .waveforms import Waveform


@dataclass(frozen=True)
class Coil:
    radius: float
    z: float
    waveform: Waveform
    turns: float = 1.0
    polarity: int = 1
    name: str = "coil"

    def __post_init__(self):
        if not all(np.isfinite(v) for v in (self.radius, self.z, self.turns)) or self.radius <= 0 or self.turns <= 0:
            raise ValueError("Coil radius/turns must be positive and geometry finite")
        if self.polarity not in (-1, 1):
            raise ValueError("Polarity must be +1 or -1")

    def current(self, t):
        return self.turns*self.polarity*self.waveform.value(t)

    def rate(self, t):
        return self.turns*self.polarity*self.waveform.derivative(t)
