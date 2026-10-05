"""Sparse exact uniform-flow linearization of the regular conservative fluid."""
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix, eye, kron


def periodic_derivative(points, spacing, coefficients):
    """Negative centered periodic first derivative with the declared stencil."""
    rows, columns, values = [], [], []
    for j in range(points):
        for offset, c in enumerate(coefficients, 1):
            rows.extend((j, j)); columns.extend(((j+offset)%points, (j-offset)%points))
            values.extend((-c/spacing, c/spacing))
    return coo_matrix((values, (rows, columns)), shape=(points, points)).tocsr()


def radial_matrices(d):
    shape = (d.elements, d.degree+1); size = np.prod(shape)
    x = d.x[..., 0]; weight = d.radial_weight[..., 0]; sign = np.array([-1., 1.])
    mass, pressure = [], []
    for column in np.eye(size):
        p = column.reshape(shape); q = p.copy(); q[0, 0] = d.sbp.axis_trace@q[0, 1:]
        nt = -2*np.einsum('ij,ej->ei', d.sbp.q, x*q)
        qt = -np.einsum('ij,eij->ei', d.sbp.q, p[:, None, :]-p[:, :, None])
        flux = np.zeros_like(q[:, [0, -1]]); star = p[:, [0, -1]].copy()
        flux[:-1, 1] = flux[1:, 0] = x[:-1, -1]*(q[:-1, -1]+q[1:, 0])
        star[:-1, 1] = star[1:, 0] = .5*(p[:-1, -1]+p[1:, 0])
        nt[:, [0, -1]] -= sign*(flux-2*x[:, [0, -1]]*q[:, [0, -1]])
        qt[:, [0, -1]] -= sign*(star-p[:, [0, -1]])
        nt /= weight; qt /= weight; qt[0, 0] = 0.
        mass.append(nt.ravel()); pressure.append(qt.ravel())
    return csr_matrix(np.column_stack(mass)), csr_matrix(np.column_stack(pressure))


def operator(fluid, background):
    d = fluid; f = d.fields(background); n = f['n'].flat[0]
    s = f['entropy'][0, 0, 0]; h = f['specific_internal'][0, 0, 0]
    v = f['velocity'][0, 0, 0, 2]
    if not (np.all(f['n'] == n) and np.all(f['entropy'] == s)
            and np.all(f['velocity'][..., :2] == 0) and np.all(f['velocity'][..., 2] == v)):
        raise ValueError('Uniform physical axial-flow background required')
    pressure = np.zeros(6)
    pressure[0] = np.sum((d.thermal.gammas-1)*h*(d.thermal.gammas-s))
    pressure[4:] = (d.thermal.gammas-1)*h
    kr = np.zeros((6, 6)); kr[[0, 3, 4, 5], 1] = [1., v, *s]
    kp = np.zeros((6, 6)); kp[1] = pressure/d.mass
    kz = np.zeros((6, 6)); kz[0, 3] = 1.; kz[1, 1] = kz[2, 2] = v
    kz[3, 3] = 2*v; kz[3, 0] = -v*v; kz[3] += pressure/d.mass
    for k in range(2):
        kz[4+k, 0] = -s[k]*v; kz[4+k, 3] = s[k]; kz[4+k, 4+k] = v
    mass, gradient = radial_matrices(d)
    z = periodic_derivative(d.nz, d.dz, d.coefficients)
    return (kron(kron(mass, eye(d.nz)), kr)+kron(kron(gradient, eye(d.nz)), kp)
            +kron(kron(eye(mass.shape[0]), z), kz)).tocsr()
