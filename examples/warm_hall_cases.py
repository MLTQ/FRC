"""Declared local energy, boundary and coupled warm Hall verification cases."""
from dataclasses import asdict, replace
from scipy.constants import atomic_mass, elementary_charge as E
from .electron_reference import ElectronReference
from .warm_hall_reference import WarmHallWave


def protocol():
    cases = []
    thermal = ElectronReference()
    wave = WarmHallWave(.035, .24, 1e20, 39.948*atomic_mass, .005, 6e20*E, amplitude=1.)
    def add(name, kind, nr, ref=None, **kw):
        cases.append(dict(name=name, kind=kind, nr=nr, nz=2*nr, reference=asdict(ref) if ref else {},
                          tolerance=1e-10, quadrature=5, steps=None, hall_iterations=16, pressure_iterations=80,
                          fluid_cfl=.25, fluid_order=2, positivity_reconstruction=True, **kw))
    for nr in (8, 16, 32, 64, 128):
        add(f'pressure_spatial_{nr}', 'pressure_spatial', nr, thermal, barotropic=False)
        add(f'trace_diffusion_{nr}', 'trace_diffusion', nr, thermal, resistivity=2e-5, decay_exponent=.02)
        cases[-1]['nz'] = 4
    add('trace_diffusion_half_time', 'trace_diffusion', 128, thermal, resistivity=2e-5, decay_exponent=.02)
    cases[-1].update(nz=4, steps=2)
    for label, ref, barotropic in (('gamma14', replace(thermal, gamma=1.4), False),
                                   ('uniform', replace(thermal, density_contrast=0.), False),
                                   ('barotropic', thermal, True)):
        for nr in (16, 32, 64):
            add(f'pressure_{label}_{nr}', 'pressure_spatial', nr, ref, barotropic=barotropic)
    for nr in (8, 16, 32, 64):
        add(f'pressure_evolve_{nr}', 'pressure_evolve', nr, thermal, end_time=1e-6, ion_pressure=16.)
        cases[-1]['steps'] = 80
    for count in (10, 20, 40, 320):
        add(f'pressure_time_{count}', 'pressure_evolve', 8, thermal, end_time=1e-6, ion_pressure=16.)
        cases[-1]['steps'] = count
    for branch in (0, 1, 2):
        ref = replace(wave, branch=branch)
        for nr in (4, 8, 16, 32):
            add(f'wave_b{branch}_{nr}', 'wave', nr, ref, phase=.5, phase_step=.025)
        add(f'wave_reverse_b{branch}', 'wave', 16, replace(ref, axial_mode=-1), phase=.5, phase_step=.025)
    for count in (8, 16, 32, 256):
        add(f'wave_time_{count}', 'wave', 8, wave, phase=.5, phase_step=.5/count)
        cases[-1]['steps'] = count
    for amp in (.1, 10., 100.):
        add(f'wave_amplitude_{amp:g}', 'wave', 16, replace(wave, amplitude=amp), phase=.5, phase_step=.025)
    add('wave_cold_limit', 'wave', 16, replace(wave, electron_fraction=1e-6), phase=.5, phase_step=.025)
    initial = dict(ion_mass=39.948*atomic_mass, bias=.005, poloidal_amplitude=.002,
                   wall_density_contrast=.1, electron_pressure=60., pressure_contrast=.2,
                   wall_pressure_contrast=.1, ion_pressure=16.)
    for nr in (8, 16, 32):
        add(f'nonlinear_{nr}', 'nonlinear', nr, thermal, initial=initial, end_time=1e-6, dt_max=2.5e-8)
    for count in (10, 20, 40, 160):
        add(f'nonlinear_time_{count}', 'nonlinear', 8, thermal, initial=initial, end_time=1e-6, dt_max=1e-6/count)
        cases[-1]['steps'] = count
    return cases
