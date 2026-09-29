"""Declared compatible Hall spatial and coupled-evolution verification cases."""
from dataclasses import asdict, replace
from scipy.constants import atomic_mass
from .hall_manufactured import HallManufactured
from .compressible_hall_reference import CompressibleHallWave


def protocol():
    cases = []
    nonlinear = HallManufactured()
    wave = CompressibleHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 16., amplitude=1.)
    def add(name, kind, ref, nr, **extra):
        config = dict(name=name, kind=kind, reference=asdict(ref), nr=nr, nz=2*nr,
                      quadrature=5, resistivity=0., nonlinear_tolerance=1e-10, dt_max=None)
        cases.append(config | extra)
    for nr in (8, 16, 32, 64, 128):
        add(f'spatial_n{nr}', 'spatial', nonlinear, nr)
    for label, changes in (('uniform', dict(density_contrast=0.)), ('strong', dict(density_contrast=.8)),
                            ('no_toroidal', dict(toroidal_slope=0.)), ('no_poloidal', dict(poloidal_amplitude=0.))):
        for nr in (16, 32, 64):
            add(f'spatial_{label}_n{nr}', 'spatial', replace(nonlinear, **changes), nr)
    for order in (3, 7):
        add(f'spatial_q{order}', 'spatial', nonlinear, 32, quadrature=order)
    for nr in (8, 16, 32, 64):
        add(f'spatial_radial_n{nr}', 'spatial', nonlinear, nr, nz=256)
    for branch in (0, 1, 2):
        phase = .5 if branch == 0 else 2*3.141592653589793
        for nr in (8, 16, 32):
            add(f'wave_b{branch}_n{nr}', 'wave', replace(wave, branch=branch), nr, phase=phase)
        add(f'wave_negative_b{branch}', 'wave', replace(wave, branch=branch, axial_mode=-1), 16, phase=phase)
    for amplitude in (10., 100.):
        add(f'wave_amplitude_{int(amplitude)}', 'wave', replace(wave, amplitude=amplitude), 16, phase=2*3.141592653589793)
    add('wave_mhd_control', 'wave', replace(wave, hall=False), 16, phase=2*3.141592653589793)
    for count in (8, 16, 32, 256):
        add(f'wave_time_{count}', 'wave', wave, 8, phase=1., dt_max=1/(abs(wave.omega)*count))
    for tolerance in (1e-8, 1e-12):
        add(f'wave_tolerance_{tolerance}', 'wave', wave, 16, phase=2*3.141592653589793, nonlinear_tolerance=tolerance)
    for nr in (8, 16, 32, 64):
        add(f'nonlinear_n{nr}', 'nonlinear', nonlinear, nr, end_time=1e-6, dt_max=1e-7)
    add('nonlinear_resistive', 'nonlinear', nonlinear, 16, end_time=1e-6, dt_max=1e-7, resistivity=1e-4)
    add('nonlinear_strong_density', 'nonlinear', replace(nonlinear, density_contrast=.8), 16, end_time=1e-6, dt_max=1e-7)
    for count in (4, 8, 16, 128):
        add(f'nonlinear_time_{count}', 'nonlinear', nonlinear, 8, end_time=3e-7, dt_max=3e-7/count)
    return cases
