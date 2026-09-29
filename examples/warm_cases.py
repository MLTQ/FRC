"""Declared spatial, time, amplitude and thermodynamic controls for warm flow."""
from dataclasses import asdict
from .warm_reference import WarmReference


def protocol():
    records = []
    def add(name, kind, nz, **changes):
        reference = asdict(WarmReference(kind=kind))
        reference.update(changes.pop('reference', {}))
        case = dict(name=name, nr=4, nz=nz, radius=.035, reference=reference,
                    periods=1. if kind == 'contact' else .5, cfl=.25, order=2,
                    density_scale=1e20, fixed_steps=None)
        case.update(changes)
        records.append(case)
    for kind in ('contact', 'simple_wave'):
        for nz in (16, 32, 64, 128, 256):
            add(f'{kind}_{nz}', kind, nz)
        add(f'{kind}_half_cfl_128', kind, 128, cfl=.125)
        add(f'{kind}_gamma14_64', kind, 64, reference=dict(gamma=1.4))
        add(f'{kind}_first_order_128', kind, 128, order=1)
        for steps in (128, 256, 512, 2048):
            add(f'{kind}_time_{steps}', kind, 32, fixed_steps=steps, periods=.05)
    add('contact_reverse_64', 'contact', 64, reference=dict(velocity=-2400.))
    add('contact_strong_64', 'contact', 64, reference=dict(density_amplitude=.6, temperature_amplitude=.16))
    add('contact_scale_64', 'contact', 64, density_scale=7e18)
    add('contact_radial_8', 'contact', 64, nr=8)
    add('contact_cold_limit_64', 'contact', 64, reference=dict(electron_ev=1e-6, ion_ev=5.999999))
    add('simple_wave_strong_64', 'simple_wave', 64, reference=dict(wave_mach=.12))
    return records
