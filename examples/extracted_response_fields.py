"""Independent physical tangent and localization of a saved response."""
import numpy as np
from plasma_bench.regular_current import RegularCurrent


def fields(sample, system, base, du, db):
    n, _, _, dn, velocity = sample.fluid(base, du)
    entropy0 = np.stack([sample.scalar(base[..., k])/n for k in (4, 5)], axis=-1)
    entropy = np.stack([(sample.scalar(du[..., k])-entropy0[..., k-4]*dn)/n for k in (4, 5)], axis=-1)
    model = RegularCurrent(system.mapping, base[..., 0], system.ion_mass)
    return dict(density=dn[..., None], velocity=velocity, entropy=entropy,
                magnetic=sample.magnetic(db), current=sample.electric(model.current(db)))


def inner(sample, a, b):
    return float(np.sum(sample.weight[..., None]*a*b))


def localization(sample, value):
    """Physical marginal power and Fourier integrals, with uncaptured norm."""
    radial_weight = sample.weight.sum(axis=1)/sample.d.length
    axial_weight = sample.weight.sum(axis=0)/(np.pi*sample.d.radius**2)
    total = inner(sample, value, value)
    radial_power = np.sum(axial_weight[None, :, None]*value**2, axis=(1, 2))
    wave = np.arange(-2*sample.mesh.nz, 2*sample.mesh.nz+1)
    phase = np.exp(-2j*np.pi*sample.z[:, None]*wave[None, :]/sample.d.length)
    coefficients = np.einsum('rzc,zk->rkc', value,
        axial_weight[:, None]*phase/sample.d.length, optimize=True)
    spectral_power = sample.d.length*np.sum(radial_weight[:, None, None]*abs(coefficients)**2, axis=(0, 2))
    summary = dict(component_power_fractions=(np.sum(sample.weight[..., None]*value**2, axis=(0, 1))/total).tolist(),
        first_element_power_fraction=float(np.sum(radial_weight*(sample.r < sample.d.radius/sample.d.elements)*radial_power)/total),
        last_element_power_fraction=float(np.sum(radial_weight*(sample.r > sample.d.radius*(1-1/sample.d.elements))*radial_power)/total),
        outer_field_cell_power_fraction=float(np.sum(radial_weight*(sample.r > sample.mesh.rf[-2])*radial_power)/total),
        maximum_fourier_index=int(wave[np.argmax(spectral_power)]),
        above_quarter_grid_fourier_power_fraction=float(spectral_power[abs(wave)>sample.mesh.nz/4].sum()/total),
        fourier_uncaptured_power_fraction=float(1-spectral_power.sum()/total))
    return summary, dict(radial_power=radial_power, fourier_index=wave, fourier_power=spectral_power)
