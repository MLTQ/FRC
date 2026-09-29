"""Declared electron-pressure spatial, drift, battery and time controls."""
from dataclasses import asdict, replace
from .electron_reference import ElectronReference


def protocol():
    cases = []
    ref = ElectronReference()
    def add(name, kind, nr, **options):
        cases.append(dict(name=name, kind=kind, nr=nr, nz=2*nr, reference=asdict(ref),
                          steps=40*max(1, nr//32), end_time=1e-6, tolerance=1e-10,
                          density_scale=1e20, expected_failure=False) | options)
    for nr in (8, 16, 32, 64, 128):
        add(f'spatial_n{nr}', 'spatial', nr)
    for label, changes in [('uniform', dict(density_contrast=0.)), ('isothermal', dict(temperature_contrast=0.)),
                            ('strong', dict(density_contrast=.8)), ('gamma', dict(gamma=1.4))]:
        for nr in (16, 32, 64):
            add(f'spatial_{label}_n{nr}', 'spatial', nr, reference=asdict(replace(ref, **changes)))
    for nz in (16, 32, 64, 128):
        add(f'drift_nz{nz}', 'drift', 4, nz=nz, steps=4*nz, amplitude=.3, slope=.01)
    add('drift_strong', 'drift', 4, nz=64, steps=256, amplitude=.7, slope=.01)
    add('drift_negative', 'drift', 4, nz=32, steps=128, amplitude=.3, slope=-.01)
    for steps in (64, 128, 256, 1024):
        add(f'drift_time_{steps}', 'drift', 4, nz=32, steps=steps, amplitude=.3, slope=.01)
    for kind in ('battery', 'nonlinear'):
        for nr in (8, 16, 32, 64):
            add(f'{kind}_n{nr}', kind, nr)
    for steps in (10, 20, 40, 320):
        add(f'nonlinear_time_{steps}', 'nonlinear', 8, steps=steps)
    for nr in (16, 32, 64):
        add(f'battery_strong_n{nr}', 'battery', nr, reference=asdict(replace(ref, density_contrast=.8, temperature_contrast=.6)), steps=160)
    add('nonlinear_scale', 'nonlinear', 16, density_scale=7e18)
    for tolerance in (1e-8, 1e-12):
        add(f'nonlinear_tolerance_{tolerance}', 'nonlinear', 16, tolerance=tolerance)
    add('battery_roundoff_control', 'battery', 64, tolerance=2e-13, expected_failure=True)
    return cases
