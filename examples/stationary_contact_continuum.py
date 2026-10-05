"""Separate Fourier/Bessel linearization of original two-species equations."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.special import jn_zeros, j0, j1
from .finite_mass_contact_reference import FiniteMassContact
from .finite_mass_wave_reference import FiniteMassWave


def operator(points=32, radial_mode=1, profile='contact'):
    ref = FiniteMassContact(velocity=0.)
    count = points; z = np.arange(count)*ref.length/count-.5*ref.length
    k = 2*np.pi*np.fft.fftfreq(count, ref.length/count); k[count//2] = 0.
    derivative = np.fft.ifft(1j*k[:, None]*np.fft.fft(np.eye(count), axis=0), axis=0).real
    q = float(jn_zeros(1, radial_mode)[-1]/ref.radius)
    f = ref.fields(z); n = f['number_density']; pressures = [f['electron_pressure'], f['ion_pressure']]
    fundamental = 2*np.pi/ref.length
    dn = -ref.number_density*ref.density_amplitude*fundamental*np.sin(fundamental*z)
    dpe = ref.number_density*E*ref.electron_ev*ref.pressure_amplitude*fundamental*np.cos(fundamental*z+.23)
    dp = [dpe, -dpe]
    if profile == 'uniform':
        n = np.full(count, ref.number_density)
        pressures = [np.full(count, ref.number_density*E*t) for t in (ref.electron_ev, ref.ion_ev)]
        dn = np.zeros(count); dp = [dn, dn]
    elif profile != 'contact': raise ValueError('Known original stationary profile required')
    mass = ref.ion_mass+ref.electron_mass
    alpha = ref.electron_mass/mass; beta = ref.ion_mass/mass
    reduced_mass = ref.ion_mass*ref.electron_mass/mass
    identity = np.eye(count); zero = np.zeros_like(identity)
    curl = np.block([[zero, -derivative, zero], [derivative, zero, q*identity], [zero, q*identity, zero]])
    current = curl/MU0
    pressure = []; entropy_gradient = []
    for species, gamma in enumerate((ref.gamma_e, ref.gamma_i)):
        block = np.zeros((count, 9*count))
        block[:, :count] = np.diag(gamma*pressures[species]/n)
        block[:, (4+species)*count:(5+species)*count] = np.diag(pressures[species])
        pressure.append(block)
        entropy_gradient.append(dp[species]/pressures[species]-gamma*dn/n)
    total_pressure = sum(pressure)
    matrix = np.zeros((9*count, 9*count)); columns = slice(6*count, 9*count)
    matrix[:count, count:2*count] = -q*np.diag(n)
    matrix[:count, 3*count:4*count] = -derivative@np.diag(n)
    matrix[count:2*count] = np.diag(q/(mass*n))@total_pressure
    matrix[count:2*count, columns] += np.diag(ref.bias/(mass*n))@current[count:2*count]
    matrix[2*count:3*count, columns] = -np.diag(ref.bias/(mass*n))@current[:count]
    matrix[3*count:4*count] = -np.diag(1/(mass*n))@derivative@total_pressure
    pressure_curl = np.zeros((3*count, 9*count))
    for species, fraction in enumerate((-beta, alpha)):
        rows = slice((4+species)*count, (5+species)*count)
        matrix[rows, 3*count:4*count] = -np.diag(entropy_gradient[species])
        matrix[rows, columns] = -np.diag(fraction*entropy_gradient[species]/(E*n))@current[2*count:]
        pressure_curl[count:2*count] -= np.diag(fraction*q*dn/(E*n*n))@pressure[species]
        pressure_curl[count:2*count, :count] += np.diag(fraction*q*dp[species]/(E*n*n))
    forcing = np.zeros((3*count, 9*count))
    forcing[:count, 2*count:3*count] = -ref.bias*identity
    forcing[count:2*count, count:2*count] = ref.bias*identity
    forcing[:count, columns] = np.diag((beta-alpha)*ref.bias/(E*n))@current[count:2*count]
    forcing[count:2*count, columns] = -np.diag((beta-alpha)*ref.bias/(E*n))@current[:count]
    inertia = np.eye(3*count)+curl@np.diag(np.tile(reduced_mass/(E*E*n), 3))@current
    matrix[6*count:] = np.linalg.solve(inertia, -curl@forcing+pressure_curl)
    speed = np.sqrt((ref.gamma_e*ref.electron_ev+ref.gamma_i*ref.ion_ev)*E/mass)
    magnetic_scale = np.sqrt(MU0*mass*ref.number_density)*speed
    scale = np.repeat([ref.number_density, speed, speed, speed, 1., 1., magnetic_scale, magnetic_scale, magnetic_scale], count)
    scaled = matrix*scale[None, :]/scale[:, None]
    return scaled, dict(z=z, q=np.array(q), scale=scale, density=n,
        pressures=np.array(pressures), density_gradient=dn, pressure_gradients=np.array(dp),
        entropy_gradients=np.array(entropy_gradient), derivative=derivative, curl=curl, inertia=inertia)


def uniform_wave_errors(points=32):
    matrix, arrays = operator(points, 1, 'uniform'); ref = FiniteMassContact()
    errors = []
    for branch in range(3):
        wave = FiniteMassWave(radius=ref.radius, length=ref.length, number_density=ref.number_density,
            ion_mass=ref.ion_mass, bias=ref.bias, pressure=ref.total_pressure,
            gamma=ref.gamma_e, gamma_i=ref.gamma_i, electron_fraction=ref.electron_ev/(ref.electron_ev+ref.ion_ev), branch=branch)
        radius = ref.radius/3; f = wave.fields(radius, arrays['z'], complex_output=True)
        n = f['density']/(ref.ion_mass+ref.electron_mass)/j0(wave.q*radius)
        v = (f['velocity']/np.array([j1(wave.q*radius), j1(wave.q*radius), j0(wave.q*radius)])).T
        magnetic = (f['magnetic']/np.array([j1(wave.q*radius), j1(wave.q*radius), j0(wave.q*radius)])).T
        entropy = [f[key]/(ref.number_density*E*t*j0(wave.q*radius))-gamma*n/ref.number_density
            for key, t, gamma in [('electron_pressure', ref.electron_ev, ref.gamma_e), ('ion_pressure', ref.ion_ev, ref.gamma_i)]]
        value = np.vstack((n, v, entropy, magnetic)).ravel()/arrays['scale']
        expected = -1j*wave.omega*value
        errors.append(float(np.linalg.norm(matrix@value-expected)/np.linalg.norm(expected)))
    return errors
