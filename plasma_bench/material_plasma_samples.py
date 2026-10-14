"""Exact material-basis fields at interior reference labels, for physical plots."""
import numpy as np
from scipy.constants import elementary_charge as E
from .compatible_polynomials import bases


def material_samples(model, state, reference_r, reference_z, reference_density, *, time=0., response=None):
    """Sample saved material fields; the caller supplies its initial density law."""
    r, z, n0 = np.broadcast_arrays(np.asarray(reference_r, float), np.asarray(reference_z, float),
                                   np.asarray(reference_density, float))
    shape = r.shape; r = r.ravel(); z = z.ravel(); n0 = n0.ravel()
    mesh = model.space.mesh; fluid = model.fluid
    if not np.isfinite(r+z+n0).all() or np.any(n0 <= 0) or np.any((r <= 0) | (r >= mesh.radius) | (abs(z) >= mesh.half_length)):
        raise ValueError('Positive density and strictly interior material labels required; no extrapolation')
    i = np.searchsorted(mesh.rf, r, side='right')-1
    j = np.searchsorted(mesh.zf, z, side='right')-1
    u = (r-mesh.rf[i])/mesh.dr; v = (z-mesh.zf[j])/mesh.dz
    weights = np.column_stack(((1-u)*(1-v), u*(1-v), u*v, (1-u)*v))
    derivative = np.stack((np.column_stack((-(1-v), 1-v, v, -v))/mesh.dr,
                           np.column_stack((-(1-u), -u, u, 1-u))/mesh.dz), axis=-1)
    ids = fluid.ids[i*mesh.nz+j]; coordinates = state.fluid.coordinates[ids]
    points = np.einsum('pn,pna->pa', weights, coordinates)
    poloidal = np.einsum('pnb,pna->pab', derivative, coordinates)
    twist = np.einsum('pna,pn->pa', derivative, state.fluid.angle[ids])
    deformation = np.zeros((len(r), 3, 3))
    deformation[:, 0, 0] = poloidal[:, 0, 0]; deformation[:, 0, 2] = poloidal[:, 0, 1]
    deformation[:, 2, 0] = poloidal[:, 1, 0]; deformation[:, 2, 2] = poloidal[:, 1, 1]
    deformation[:, 1, 0] = points[:, 0]*twist[:, 0]
    deformation[:, 1, 1] = points[:, 0]/r
    deformation[:, 1, 2] = points[:, 0]*twist[:, 1]
    determinant = np.linalg.det(deformation)
    if np.any(determinant <= 0):
        raise ValueError('Nonpositive actual material geometry')
    nx, dx, ex, _ = bases(.5*mesh.rf**2, .5*r*r, model.space.degree)
    nz, dz, ez, _ = bases(mesh.zf, z, model.space.degree)
    nx, dx, ex, nz, dz, ez = [a.toarray() for a in (nx, dx, ex, nz, dz, ez)]
    pair = lambda left, coefficient, right: np.sum((left@coefficient)*right, axis=1)
    nt = model.space.toroidal_magnetic_count
    toroidal = state.magnetic[:nt].reshape(mesh.nr, mesh.nz)
    poloidal_flux = state.magnetic[nt:].reshape(mesh.nr, mesh.nz+1)
    flux = pair(nx[:, 1:], poloidal_flux, nz)
    b0 = np.column_stack((-pair(nx[:, 1:], poloidal_flux, dz)/r,
                           r*pair(ex, toroidal, ez), pair(dx[:, 1:], poloidal_flux, nz)))
    induced = np.einsum('pij,pj->pi', deformation, b0)/determinant[:, None]
    if response is None:
        response = model.response(state, rate=False, time=time) if hasattr(model, 'source') else model.response(state, rate=False)
    w = response['drift_coefficients']; er_count = mesh.nr*(mesh.nz+1)
    ez_end = model.space.meridional_electric_count
    wr = w[:er_count].reshape(mesh.nr, mesh.nz+1)
    wz = w[er_count:ez_end].reshape(mesh.nr+1, mesh.nz)
    wt = w[ez_end:].reshape(mesh.nr, mesh.nz+1)
    w0 = np.column_stack((r*pair(ex, wr, nz), pair(nx[:, 1:], wt, nz)/r, pair(nx, wz, ez)))
    drift = np.linalg.solve(deformation.transpose(0, 2, 1), w0[..., None])[..., 0]
    density = n0/determinant; current = E*density[:, None]*drift
    magnetic = induced.copy()
    if hasattr(model, 'source'):
        potential, unit = model.source(state.fluid.coordinates, points)
        amplitude = float(model.waveform.value(time))
        magnetic += amplitude*unit; flux += amplitude*points[:, 0]*potential[:, 1]
    fields = response['fields']
    meridional = np.einsum('pn,pna->pa', weights, fields['velocity'][ids])
    omega = np.einsum('pn,pn->p', weights, fields['omega'][ids])
    velocity = np.column_stack((meridional[:, 0], points[:, 0]*omega, meridional[:, 1]))
    return dict(points=points.reshape(shape+(2,)), magnetic=magnetic.reshape(shape+(3,)),
        induced_magnetic=induced.reshape(shape+(3,)), current=current.reshape(shape+(3,)),
        velocity=velocity.reshape(shape+(3,)), flux=flux.reshape(shape),
        density=density.reshape(shape), volume_ratio=determinant.reshape(shape))
