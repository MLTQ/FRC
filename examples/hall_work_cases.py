"""Declared spatial Hall work and finite-pressure mode verification protocol."""
from dataclasses import asdict, replace
from scipy.constants import atomic_mass
from .hall_manufactured import HallManufactured
from .compressible_hall_reference import CompressibleHallWave


def protocol():
    cases = []
    base = HallManufactured()
    def add(name, kind, reference, nr, quadrature=3, density='quadrature'):
        cases.append(dict(name=name, kind=kind, reference=asdict(reference), nr=nr, nz=2*nr,
                          quadrature=quadrature, density_sampling=density))
    for nr in (8, 16, 32, 64, 128):
        add(f'nonlinear_n{nr}', 'nonlinear', base, nr)
    for name, change in (('uniform', dict(density_contrast=0.)), ('strong_density', dict(density_contrast=.8)),
                          ('zero_toroidal', dict(toroidal_slope=0.)), ('zero_poloidal', dict(poloidal_amplitude=0.))):
        for nr in (16, 32, 64):
            add(f'{name}_n{nr}', 'nonlinear', replace(base, **change), nr)
    for nr, order in ((32, 5), (32, 7), (64, 5)):
        add(f'quadrature_q{order}_n{nr}', 'nonlinear', base, nr, order)
    for nr in (32, 64):
        add(f'cell_density_n{nr}', 'nonlinear', base, nr, density='cell')
    wave = CompressibleHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 16., amplitude=.01)
    for branch in (0, 1, 2):
        for sign in (-1, 1):
            for nr in (8, 16, 32):
                add(f'wave_b{branch}_k{sign}_n{nr}', 'wave', replace(wave, branch=branch, axial_mode=sign), nr)
    for amplitude in (.001, .1):
        add(f'wave_amplitude_{amplitude}', 'wave', replace(wave, amplitude=amplitude), 32)
    return cases
