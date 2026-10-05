"""Matched physical Fourier/Bessel projection of a production mode."""
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.constants import atomic_mass, electron_mass, elementary_charge as E
from .stationary_contact_continuum import operator


def constant_radial_operator(points):
    """Original radially constant scalar/axial acoustic and entropy equations."""
    _, data = operator(points, 1)
    n = data['density']; pe, pi = data['pressures']; derivative = data['derivative']
    mass = 39.948*atomic_mass+electron_mass
    pressure = np.zeros((points, 9*points))
    pressure[:, :points] = np.diag(((5/3)*pe+1.4*pi)/n)
    pressure[:, 4*points:5*points] = np.diag(pe); pressure[:, 5*points:6*points] = np.diag(pi)
    matrix = np.zeros((9*points, 9*points))
    matrix[:points, 3*points:4*points] = -derivative@np.diag(n)
    matrix[3*points:4*points] = -np.diag(1/(mass*n))@derivative@pressure
    for k in range(2): matrix[(4+k)*points:(5+k)*points, 3*points:4*points] = -np.diag(data['entropy_gradients'][k])
    return matrix*data['scale'][None, :]/data['scale'][:, None], data


def project(fields, points=64, modes=16):
    r, z = fields['radius'], fields['axial']; wr, wz = fields['radial_weight'], fields['axial_weight']
    physical = fields['physical_mode']; radius = .035; length = .24
    frequency = 2*np.pi*np.fft.fftfreq(points, length/points)
    transform = np.exp(-1j*frequency[:, None]*(z[None, :]+length/2))*wz[None, :]/length
    # Use complete conjugate pairs; omit the unpaired even-grid Nyquist mode.
    transform[points//2] = 0.
    target_z = np.arange(points)*length/points-length/2
    evaluate = np.exp(1j*(target_z[:, None]+length/2)*frequency[None, :])
    scalar = np.array([True, False, False, True, True, True, False, False, True])
    qs = np.r_[0., jn_zeros(1, modes)/radius]
    projected = []; predicted = []; captured = []; norms = []
    scales = operator(points, 1)[1]['scale'].reshape(9, points)[:, 0]
    original_norms = np.sum(wr[:, None, None]*wz[None, :, None]*(physical/scales)**2, axis=(0, 1))
    for radial_mode, q in enumerate(qs):
        if radial_mode:
            basis = np.where(scalar[None, :], j0(q*r)[:, None], j1(q*r)[:, None])
            norm = np.pi*radius**2*j0(q*radius)**2
            matrix, data = operator(points, radial_mode)
        else:
            basis = np.broadcast_to(scalar, (len(r), 9)).astype(float); norm = np.pi*radius**2
            matrix, data = constant_radial_operator(points)
        amplitude = np.einsum('r,rk,rzk->kz', wr, basis, physical)/norm
        values = (evaluate@transform@amplitude.T).T
        if np.max(abs(values.imag)) > 1e-10*max(np.max(abs(values.real)), np.finfo(float).tiny):
            raise ValueError('Nonreal matched Fourier projection')
        scaled = values.real.ravel()/data['scale']
        projected.append(scaled.reshape(9, points)); predicted.append((matrix@scaled).reshape(9, points))
        captured.append(norm*length*np.mean(scaled.reshape(9, points)**2, axis=1)); norms.append(norm)
    return dict(projected=np.array(projected), predicted=np.array(predicted), radial_norm=np.array(norms),
        original_norms=original_norms, captured_norms=np.sum(captured, axis=0), q=qs, axial=target_z, scales=scales)
