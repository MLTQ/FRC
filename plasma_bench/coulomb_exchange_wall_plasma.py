"""Conditional Maxwellian temperature exchange in the actual moving state."""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0
from .resistive_wall_plasma import ResistiveWallPlasma


def equilibration_rate(number_density, electron_ev, ion_ev, ion_mass):
    """Z=1, small relative drift; NRL electron-ion log, SI Landau rate."""
    n, te, ti = np.broadcast_arrays(number_density, electron_ev, ion_ev)
    if (not all(np.isfinite(x).all() for x in (n, te, ti))
            or any(np.any(x <= 0) for x in (n, te, ti))
            or not np.isfinite(ion_mass) or ion_mass <= ME):
        raise ValueError('Positive plasma state and physical ion mass required')
    if np.any(ti*ME/ion_mass >= np.minimum(te, 10.)):
        raise ValueError('Unsupported electron-ion Coulomb logarithm branch')
    log = np.where(te < 10., 23.-.5*np.log(n/1e6)+1.5*np.log(te),
        24.-.5*np.log(n/1e6)+np.log(te))
    if np.any(log <= 0):
        raise ValueError('Nonpositive Coulomb logarithm; no floor applied')
    coupling = E**4/(4*np.pi*EPS0)**2
    thermal_speed2 = E*(te/ME+ti/ion_mass)
    rate = 8*np.sqrt(2*np.pi)*n*coupling*log/(3*ME*ion_mass*thermal_speed2**1.5)
    if not np.isfinite(rate).all() or np.any(rate <= 0):
        raise ValueError('Collision rate exceeds positive finite range')
    return rate, log


class CoulombExchangeWallPlasma(ResistiveWallPlasma):
    def __init__(self, mapping, ion_mass, **options):
        if not np.allclose(mapping.fluid.thermal.gammas, 5/3, rtol=0., atol=4*np.finfo(float).eps):
            raise ValueError('Monatomic electron and ion gamma=5/3 required for this closure')
        super().__init__(mapping, ion_mass, **options)

    def derivative(self, state, field):
        result = super().derivative(state, field)
        fields = self.fluid.fields(state)
        te, ti = np.moveaxis(fields['temperature'], -1, 0)
        frequency, log = equilibration_rate(fields['n'], te, ti, self.ion_mass)
        q = 1.5*E*fields['n']*frequency*(ti-te)
        exchange = np.stack((q, -q), axis=-1)/fields['specific_internal']
        fluid = result.fluid.copy()
        fluid[..., 4:] += exchange
        if not np.isfinite(fluid).all():
            raise ValueError('Nonfinite collisional thermal rate')
        net = float(np.sum(self.fluid.volume[..., None]*fields['specific_internal']*exchange))
        powers = dict(result.diagnostics['powers'])
        powers['fluid'] += net
        diagnostics = dict(result.diagnostics, powers=powers, total_power=sum(powers.values()),
            exchange_power_to_electrons=float(np.sum(self.fluid.volume*q)),
            exchange_net_grid_power=net,
            exchange_entropy_rate_over_kB=float(np.sum(self.fluid.volume[..., None]*exchange/(5/3-1))),
            exchange_frequency_range_s=[float(frequency.min()), float(frequency.max())],
            exchange_coulomb_log_range=[float(log.min()), float(log.max())],
            exchange_small_log_review=bool(np.any(log < 5.)))
        return replace(result, fluid=fluid, diagnostics=diagnostics)
