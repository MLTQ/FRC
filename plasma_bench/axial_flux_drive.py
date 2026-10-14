"""Prescribed enclosed axial flux and its consistent inductive source work."""
from dataclasses import dataclass, replace
import numpy as np
from scipy.constants import mu_0 as MU0
from .coulomb_exchange_wall_plasma import CoulombExchangeWallPlasma
from .spitzer_wall_plasma import SpitzerWallPlasma
from .regular_current import RegularCurrent
from .variable_compatible_current import weak_vector


@dataclass(frozen=True)
class AxialFluxPulse:
    peak_field_t: float
    duration_s: float

    def __post_init__(self):
        if not np.isfinite([self.peak_field_t, self.duration_s]).all() or self.duration_s <= 0:
            raise ValueError('Finite flux-field amplitude and positive pulse duration required')

    def value(self, when):
        t = np.asarray(when, float)
        if not np.isfinite(t).all():
            raise ValueError('Finite source time required')
        x = np.clip(t/self.duration_s, 0., 1.)
        return np.where((t > 0) & (t < self.duration_s),
            self.peak_field_t*np.sin(np.pi*x)**2, 0.)

    def derivative(self, when):
        t = np.asarray(when, float)
        if not np.isfinite(t).all():
            raise ValueError('Finite source time required')
        x = np.clip(t/self.duration_s, 0., 1.)
        return np.where((t > 0) & (t < self.duration_s),
            self.peak_field_t*np.pi/self.duration_s*np.sin(2*np.pi*x), 0.)


class DrivenTransport:
    def __init__(self, transport, rate):
        self.transport = transport
        self.rate = rate

    def pair(self, state, current, density_rate):
        result = self.transport.pair(state, current, density_rate)
        geometry = self.transport.model.geometry
        external = np.zeros(geometry.r.shape+(3,))
        external[..., 1] = -.5*geometry.r*self.rate
        # E_induced=E_physical-E_external drives only the zero-wall-flux increment.
        result['forcing'] -= weak_vector(geometry, external)
        return result


class AxialFluxDrive:
    """Mixin for sequential stage evaluations of the existing moving models."""
    def __init__(self, mapping, ion_mass, *, drive, **options):
        if not callable(getattr(drive, 'value', None)) or not callable(getattr(drive, 'derivative', None)):
            raise ValueError('Explicit flux-field value and time derivative required')
        super().__init__(mapping, ion_mass, **options)
        self.drive = drive
        self._active_drive_rate = None

    def field_value(self, when):
        value = float(self.drive.value(when))
        if not np.isfinite(value):
            raise ValueError('Finite prescribed mean axial field required')
        return value

    def operators(self, model):
        if self._active_drive_rate is None:
            raise RuntimeError('Driven operators require an explicit source time')
        transport, species, induction = super().operators(model)
        return DrivenTransport(transport, self._active_drive_rate), species, induction

    def derivative(self, state, field, when):
        if field.axial_bias != self.field_value(when):
            raise ValueError('Stored mean flux and prescribed source time must agree')
        rate = float(self.drive.derivative(when))
        if not np.isfinite(rate) or self._active_drive_rate is not None:
            raise ValueError('Finite non-reentrant source evaluation required')
        self._active_drive_rate = rate
        try:
            result = super().derivative(state, field)
        finally:
            self._active_drive_rate = None
        model = RegularCurrent(self.mapping, state[..., 0], self.ion_mass)
        j = self.geometry.sample_electric(model.current(field.packed_increment(self.geometry)))
        external_power = float(np.sum(self.geometry.weight*j[..., 1]*(-.5*self.geometry.r*rate)))
        harmonic_power = float(np.sum(model.mesh.volumes)*field.axial_bias*rate/MU0)
        powers = dict(result.diagnostics['powers'])
        powers['magnetic'] += harmonic_power
        source_power = external_power+harmonic_power
        diagnostics = dict(result.diagnostics, powers=powers, total_power=sum(powers.values()),
            source_time_s=float(when), mean_axial_field_t=field.axial_bias,
            mean_axial_field_rate_t_s=rate,
            wall_loop_voltage_v=float(-np.pi*model.mesh.radius**2*rate),
            source_plasma_power_w=external_power, source_harmonic_power_w=harmonic_power,
            source_power_w=source_power, power_balance_residual_w=sum(powers.values())-source_power,
            electric_representation='Returned electric coefficients describe E_physical-E_external; add exact E_external,theta=-r*dBbar/dt/2 for physical E.',
            drive_scope='Prescribed enclosed flux pi R^2 Bbar; filled periodic wall cylinder. Not prescribed Bz at the wall, finite external coils, or a solved exterior vacuum.')
        return replace(result, diagnostics=diagnostics)


class FluxDrivenSpitzerPlasma(AxialFluxDrive, SpitzerWallPlasma):
    pass


class FluxDrivenResistivePlasma(AxialFluxDrive, CoulombExchangeWallPlasma):
    pass
