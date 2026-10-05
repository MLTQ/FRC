"""Independent Fourier/Bessel growing contact mode and its full electric field."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.special import j0, j1
from .finite_mass_contact_reference import FiniteMassContact
from .implicit_plasma_records import sha


class ContinuumContactMode:
    def __init__(self):
        root = Path('outputs/stationary_contact_audit')
        record = json.loads((root/'summary.json').read_text())['continuum'][3]
        raw = root/record['arrays_file']
        if sha(raw) != record['arrays_sha256']: raise ValueError('Frozen independent continuum mode required')
        a = dict(np.load(raw)); self.input_sha256 = sha(raw); self.ref = FiniteMassContact(velocity=0.)
        self.q = float(a['q']); self.growth = record['maximum_real']; self.points = len(a['z'])
        vector = a['leading_vector']
        if record['leading_imaginary'] != 0 or np.max(abs(vector.imag)) != 0:
            raise ValueError('Saved real growing continuum mode required')
        vector = vector.real/np.max(abs(vector.real))
        if vector[np.argmax(abs(vector))] < 0: vector = -vector
        physical = (vector*a['scale']).reshape(9, self.points)
        self.coefficients = np.fft.fft(physical, axis=1)/self.points
        self.frequency = 2*np.pi*np.fft.fftfreq(self.points, self.ref.length/self.points)
        # The source even-grid operator defines the Nyquist derivative as zero.
        self.frequency[self.points//2] = 0.
        self.derivative = a['derivative']; self.curl = a['curl']; self.physical_nodes = physical

    def amplitudes(self, z, derivative=0):
        z = np.asarray(z); unique, inverse = np.unique(z, return_inverse=True)
        phase = np.exp(1j*(unique[:, None]+self.ref.length/2)*self.frequency[None, :])
        # Nyquist value retains its cosine; only its derivative is defined zero.
        if derivative == 0:
            phase[:, self.points//2] = np.cos(np.pi*self.points*(unique+self.ref.length/2)/self.ref.length)
        values = (phase@(self.coefficients*(1j*self.frequency)**derivative).T).real
        return values[inverse].reshape(z.shape+(9,))

    def fields(self, r, z):
        r, z = np.broadcast_arrays(r, z); a = self.amplitudes(z); az = self.amplitudes(z, 1)
        scalar = j0(self.q*r); transverse = j1(self.q*r)
        basis = np.stack([scalar, transverse, transverse, scalar, scalar, scalar, transverse, transverse, scalar], axis=-1)
        physical = a*basis
        current = np.stack((-az[..., 7]*transverse,
            (az[..., 6]+self.q*a[..., 8])*transverse, self.q*a[..., 7]*scalar), axis=-1)/MU0
        ref = self.ref; f = ref.fields(z); n = f['number_density']
        pressure = np.stack((f['electron_pressure'], f['ion_pressure']), axis=-1)
        fundamental = 2*np.pi/ref.length
        nprime = -ref.number_density*ref.density_amplitude*fundamental*np.sin(fundamental*z)
        peprime = ref.number_density*E*ref.electron_ev*ref.pressure_amplitude*fundamental*np.cos(fundamental*z+.23)
        pprime = np.stack((peprime, -peprime), axis=-1)
        gammas = np.array([ref.gamma_e, ref.gamma_i])
        dp = pressure*(gammas*a[..., :1]/n[..., None]+a[..., 4:6])
        dpz = pprime*(gammas*a[..., :1]/n[..., None]+a[..., 4:6])
        dpz += pressure*(gammas*(az[..., :1]/n[..., None]-a[..., :1]*nprime[..., None]/n[..., None]**2)+az[..., 4:6])
        mass = ref.ion_mass+ref.electron_mass; alpha = ref.electron_mass/mass; beta = ref.ion_mass/mass
        reduced = ref.ion_mass*alpha; bias = np.array([0., 0., ref.bias])
        electric = -np.cross(physical[..., 1:4], bias)+(beta-alpha)*np.cross(current, bias)/(E*n[..., None])
        electric += reduced*self.growth*current/(E**2*n[..., None])
        fractions = np.array([-beta, alpha])
        electric[..., 0] += np.sum(fractions*(-self.q*dp*transverse[..., None])/n[..., None], axis=-1)/E
        electric[..., 2] += np.sum(fractions*(dpz/n[..., None]-a[..., :1]*pprime/n[..., None]**2), axis=-1)*scalar/E
        return dict(physical=physical, current=current, electric=electric,
            psi=r*transverse*a[..., 8]/self.q)
