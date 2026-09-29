"""Independent continuum checks for finite-pressure Hall integration oracles."""
from dataclasses import replace
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from numpy.testing import assert_allclose
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.linalg import eigvalsh
from examples.compressible_hall_reference import CompressibleHallWave, generator, dispersion_terms
from examples.hall_manufactured import HallManufactured


def wave(**kwargs):
    # Conditional physical Ar+ parameters; no assertion that this is a formation regime.
    values = dict(radius=.035, length=.24, number_density=1e20, ion_mass=39.948*1.66053906892e-27,
                  bias=.005, pressure=16., amplitude=.01)
    return CompressibleHallWave(**(values | kwargs))


def derivative(function, x, h):
    return (function(x-2*h)-8*function(x-h)+8*function(x+h)-function(x+2*h))/(12*h)


@pytest.mark.parametrize('branch', [0, 1, 2])
@pytest.mark.parametrize('sign', [-1, 1])
def test_wave_satisfies_independent_maxwell_and_species_equations(branch, sign):
    mode = wave(branch=branch, axial_mode=sign, bias=sign*.005, sense=-sign)
    r, z = np.meshgrid(np.linspace(.07, .93, 17)*mode.radius,
                       np.linspace(-.47, .47, 15)*mode.length, indexing='ij')
    data = mode.fields(r, z, complex_output=True)
    dr = {key: derivative(lambda x: mode.fields(x, z, complex_output=True)[key], r, 2e-5*mode.radius)
          for key in ('magnetic', 'electric', 'pressure', 'velocity')}
    dz = {key: derivative(lambda x: mode.fields(r, x, complex_output=True)[key], z, 2e-5*mode.length)
          for key in ('magnetic', 'electric', 'pressure', 'velocity')}
    def curl(key):
        return np.stack((-dz[key][..., 1], dz[key][..., 0]-dr[key][..., 2],
                         dr[key][..., 1]+data[key][..., 1]/r), axis=-1)
    assert np.linalg.norm(curl('magnetic')-MU0*data['current'])/np.linalg.norm(MU0*data['current']) < 2e-10
    assert np.linalg.norm(-curl('electric')+1j*mode.omega*data['magnetic'])/np.linalg.norm(mode.omega*data['magnetic']) < 2e-9
    gradp = np.stack((dr['pressure'], np.zeros_like(r), dz['pressure']), axis=-1)
    bias = np.array([0., 0., mode.bias])
    electron_velocity = data['velocity']-data['current']/(mode.number_density*E)
    electron_terms = (data['electric'], np.cross(electron_velocity, bias))
    assert np.linalg.norm(sum(electron_terms))/sum(np.linalg.norm(x) for x in electron_terms) < 1e-11
    ion_terms = (-1j*mode.omega*mode.density*data['velocity'],
                 -mode.number_density*E*(data['electric']+np.cross(data['velocity'], bias)), gradp)
    assert np.linalg.norm(sum(ion_terms))/sum(np.linalg.norm(x) for x in ion_terms) < 2e-10
    divergence = dr['velocity'][..., 0]+data['velocity'][..., 0]/r+dz['velocity'][..., 2]
    assert np.linalg.norm(-1j*mode.omega*data['pressure']+mode.gamma*mode.pressure*divergence)/np.linalg.norm(mode.omega*data['pressure']) < 2e-9


def test_dispersion_all_branches_and_signed_inputs():
    for q in (0., .7, 5.):
        for k in (-3., 2.):
            for cs in (.2, 1., 20.):
                for va, dh in ((1., 0.), (1., .4), (-1., -2.)):
                    matrix = generator(q, k, cs, va, dh)
                    assert_allclose(matrix+matrix.conj().T, 0., atol=0.)
                    roots = eigvalsh(1j*matrix)
                    left, right = dispersion_terms(roots, q, k, cs, va, dh)
                    scale = max(np.hypot(q, k)*cs, np.hypot(q, k)*abs(va), abs(dh*k)*np.hypot(q, k))**6
                    assert np.max(abs(left-right))/scale < 4e-14


def test_mhd_and_parallel_and_incompressible_limits():
    q, k, cs, va = 3., 2., 1.3, .8
    total2 = q*q+k*k
    combined = total2*(cs*cs+va*va)
    product = total2*k*k*cs*cs*va*va
    fast2 = .5*(combined+np.sqrt(combined**2-4*product))
    expected = np.sort([np.sqrt(fast2), np.sqrt(product/fast2), abs(k*va)])
    assert_allclose(eigvalsh(1j*generator(q, k, cs, va, 0.))[3:], expected, rtol=2e-14)
    dh = .7
    split = abs(dh)*k*k
    fast = .5*(np.hypot(split, 2*k*va)+split)
    slow = (k*va)**2/fast
    assert_allclose(eigvalsh(1j*generator(0., k, cs, va, dh))[3:], np.sort([fast, slow, abs(k*cs)]), rtol=2e-14)
    hall = abs(dh*k)*np.sqrt(total2)
    fast = .5*(np.hypot(hall, 2*k*va)+hall)
    expected = np.array([(k*va)**2/fast, fast])
    errors = []
    for sound in (30., 60., 120.):
        actual = eigvalsh(1j*generator(q, k, sound, va, dh))[3:5]
        errors.append(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


@pytest.mark.parametrize('branch', [0, 1, 2])
def test_wall_axis_and_physical_volume_wave_energy(branch):
    mode = wave(branch=branch)
    z = np.linspace(-mode.length/2, mode.length/2, 31)
    wall = mode.fields(mode.radius, z)
    assert np.max(abs(wall['velocity'][..., :2])) < 1e-15*mode.amplitude
    assert np.max(abs(wall['magnetic'][..., :2])) < 1e-15*mode.amplitude*np.sqrt(MU0*mode.density)
    assert np.max(abs(wall['electric'][..., 1:])) < 1e-14*mode.amplitude*abs(mode.bias)
    axis = mode.fields(0., z)
    assert_allclose(axis['velocity'][..., :2], 0., atol=0.)
    assert_allclose(axis['magnetic'][..., :2], 0., atol=0.)
    assert_allclose(axis['current'][..., :2], 0., atol=0.)
    x, w = leggauss(48)
    r, zz = np.meshgrid(.5*mode.radius*(x+1), .5*mode.length*x, indexing='ij')
    weight = 2*np.pi*r*w[:, None]*w[None, :]*mode.radius*mode.length/4
    for time in (0., 1.37/abs(mode.omega)):
        fields = mode.fields(r, zz, time)
        density = .5*mode.density*np.sum(fields['velocity']**2, axis=-1)
        density += .5*np.sum(fields['magnetic']**2, axis=-1)/MU0
        density += fields['pressure']**2/(2*mode.density*mode.sound_speed**2)
        assert np.sum(weight*density) == pytest.approx(mode.perturbation_energy(), rel=1e-13)
    # It is intentionally not a resistive conducting fixture: this wall trace survives.
    assert np.max(abs(wall['current'][..., 2])) > 0


def test_nonlinear_manufactured_derivatives_and_local_hall_work():
    reference = HallManufactured()
    r, z = np.meshgrid(reference.radius*np.linspace(.01, .99, 29),
                       reference.length*np.linspace(-.45, .46, 31), indexing='ij')
    data = reference.fields(r, z)
    for key in ('magnetic', 'current', 'electric'):
        for direction, x, size in (('r', r, reference.radius), ('z', z, reference.length)):
            function = (lambda x: reference.fields(x, z)[key]) if direction == 'r' else (lambda x: reference.fields(r, x)[key])
            numerical = derivative(function, x, 2e-5*size)
            exact = data[key+'_'+direction]
            assert np.linalg.norm(numerical-exact)/np.linalg.norm(exact) < 3e-10
    curl = np.stack((-data['magnetic_z'][..., 1], data['magnetic_z'][..., 0]-data['magnetic_r'][..., 2],
                     data['magnetic_r'][..., 1]+data['magnetic'][..., 1]/r), axis=-1)
    assert_allclose(curl/MU0, data['current'], rtol=3e-13, atol=1e-10)
    check_r = np.concatenate(([0.], reference.radius*np.linspace(.001, 1., 100)))[:, None]
    check = reference.fields(check_r, z[0])
    work = np.sum(check['current']*check['electric'], axis=-1)
    work_scale = np.max(np.linalg.norm(check['current'], axis=-1)*np.linalg.norm(check['electric'], axis=-1))
    assert np.max(abs(work))/work_scale < 3e-16
    magnetic_rate = np.sum(check['magnetic']*check['induction'], axis=-1)/MU0
    assert np.max(abs(magnetic_rate+check['poynting_divergence']))/np.max(abs(magnetic_rate)) < 2e-15
    assert_allclose(check['electric'][-1], 0., atol=0.)
    assert np.linalg.norm(data['induction']-replace(reference, density_contrast=0.).fields(r, z)['induction'])/np.linalg.norm(data['induction']) > .1


def test_invalid_reference_domains_are_not_regularized():
    for kwargs in (dict(number_density=0.), dict(pressure=0.), dict(gamma=1.), dict(branch=3), dict(axial_mode=0), dict(radial_mode=0)):
        with pytest.raises(ValueError):
            wave(**kwargs)
    with pytest.raises(ValueError):
        HallManufactured(density_contrast=1.)
    with pytest.raises(ValueError):
        HallManufactured().fields(-.01, 0.)
