"""Prescribed currents in amperes and analytic derivatives in amperes/second."""
from dataclasses import dataclass, asdict
from typing import Protocol, Callable
import numpy as np


class Waveform(Protocol):
    def value(self, t): ...
    def derivative(self, t): ...


def _finite(*values):
    if not all(np.isfinite(v) for v in values):
        raise ValueError("Waveform parameters must be finite")


def _positive(*values):
    _finite(*values)
    if any(v <= 0 for v in values):
        raise ValueError("Waveform times and frequencies must be positive")


@dataclass(frozen=True)
class DC:
    amplitude: float = 0.0

    def __post_init__(self):
        _finite(self.amplitude)

    def value(self, t):
        return np.zeros_like(np.asarray(t, dtype=float)) + self.amplitude

    def derivative(self, t):
        return np.zeros_like(np.asarray(t, dtype=float))


@dataclass(frozen=True)
class Gaussian:
    amplitude: float
    center: float
    width: float

    def __post_init__(self):
        _finite(self.amplitude, self.center)
        _positive(self.width)

    def value(self, t):
        x = (np.asarray(t) - self.center) / self.width
        return self.amplitude * np.exp(-0.5 * x*x)

    def derivative(self, t):
        return -(np.asarray(t)-self.center) / self.width**2 * self.value(t)


@dataclass(frozen=True)
class Unipolar:
    amplitude: float
    start: float = 0.0
    rise: float = 5e-6
    hold: float = 10e-6
    fall: float = 15e-6

    def __post_init__(self):
        _finite(self.amplitude, self.start, self.hold)
        _positive(self.rise, self.fall)
        if self.hold < 0:
            raise ValueError("Hold duration cannot be negative")

    def _parts(self, t):
        x = np.asarray(t, dtype=float)-self.start
        a = np.clip(x/self.rise, 0, 1)
        b = np.clip((x-self.rise-self.hold)/self.fall, 0, 1)
        return a, b

    def value(self, t):
        a, b = self._parts(t)
        return self.amplitude * (0.5-0.5*np.cos(np.pi*a)) * (0.5+0.5*np.cos(np.pi*b))

    def derivative(self, t):
        a, b = self._parts(t)
        da = np.where((a > 0) & (a < 1), np.pi*np.sin(np.pi*a)/(2*self.rise), 0)
        db = np.where((b > 0) & (b < 1), -np.pi*np.sin(np.pi*b)/(2*self.fall), 0)
        return self.amplitude * (da*(0.5+0.5*np.cos(np.pi*b)) + db*(0.5-0.5*np.cos(np.pi*a)))


@dataclass(frozen=True)
class DampedSine:
    amplitude: float
    frequency: float = 5e4
    decay: float = 20e-6
    start: float = 0.0

    def __post_init__(self):
        _finite(self.amplitude, self.start)
        _positive(self.frequency, self.decay)

    def value(self, t):
        x = np.maximum(np.asarray(t)-self.start, 0)
        return self.amplitude*np.exp(-x/self.decay)*np.sin(2*np.pi*self.frequency*x)

    def derivative(self, t):
        x = np.maximum(np.asarray(t)-self.start, 0)
        w = 2*np.pi*self.frequency
        return np.where(np.asarray(t) >= self.start,
                        self.amplitude*np.exp(-x/self.decay)*(w*np.cos(w*x)-np.sin(w*x)/self.decay), 0)


@dataclass(frozen=True)
class SineBurst:
    amplitude: float
    frequency: float = 5e4
    cycles: int = 3
    start: float = 0.0

    def __post_init__(self):
        _finite(self.amplitude, self.start)
        _positive(self.frequency)
        if not isinstance(self.cycles, int) or self.cycles < 1:
            raise ValueError("Burst cycles must be a positive integer")

    def value(self, t):
        duration = self.cycles/self.frequency
        x = np.clip(np.asarray(t)-self.start, 0, duration)
        return self.amplitude*np.sin(np.pi*x/duration)**2*np.sin(2*np.pi*self.frequency*x)

    def derivative(self, t):
        duration = self.cycles/self.frequency
        x = np.clip(np.asarray(t)-self.start, 0, duration)
        w = 2*np.pi*self.frequency
        d = self.amplitude * (np.pi/duration*np.sin(2*np.pi*x/duration)*np.sin(w*x)
                              + np.sin(np.pi*x/duration)**2*w*np.cos(w*x))
        return np.where((x > 0) & (x < duration), d, 0)


@dataclass(frozen=True)
class Sampled:
    times: tuple[float, ...]
    currents: tuple[float, ...]

    def __post_init__(self):
        object.__setattr__(self, "times", tuple(self.times))
        object.__setattr__(self, "currents", tuple(self.currents))
        if len(self.times) < 2 or len(self.times) != len(self.currents):
            raise ValueError("Provide matching arrays with at least two samples")
        if not np.all(np.isfinite(self.times)) or not np.all(np.isfinite(self.currents)) or np.any(np.diff(self.times) <= 0):
            raise ValueError("Samples must be finite with strictly increasing times")

    def value(self, t):
        return np.interp(t, self.times, self.currents)

    def derivative(self, t):
        t = np.asarray(t)
        slopes = np.diff(self.currents)/np.diff(self.times)
        i = np.clip(np.searchsorted(self.times, t, side="right")-1, 0, len(slopes)-1)
        return np.where((t >= self.times[0]) & (t < self.times[-1]), slopes[i], 0)


@dataclass(frozen=True)
class Callback:
    function: Callable
    rate: Callable
    description: str

    def value(self, t):
        return np.asarray(self.function(t))

    def derivative(self, t):
        return np.asarray(self.rate(t))


_TYPES = {c.__name__: c for c in (DC, Gaussian, Unipolar, DampedSine, SineBurst, Sampled)}


def to_dict(wave):
    if type(wave).__name__ not in _TYPES:
        raise ValueError("Callbacks are script-only; sample them before reproducible export")
    return {"kind": type(wave).__name__, **asdict(wave)}


def from_dict(data):
    data = dict(data)
    kind = data.pop("kind")
    if kind not in _TYPES:
        raise ValueError(f"Unknown waveform {kind!r}")
    return _TYPES[kind](**data)
