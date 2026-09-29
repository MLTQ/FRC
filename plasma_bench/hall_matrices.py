"""Sparse compatible magnetic Hessian and local current-to-electric transport."""
import numpy as np
from scipy.sparse import coo_matrix, diags, block_diag, bmat
from scipy.constants import elementary_charge as E, mu_0 as MU0
from .compatible_magnetic import radial_mass


def magnetic_hessian(mesh):
    """Exact Hessian of the fixed compatible magnetic energy on periodic DOFs."""
    nr, nz = mesh.nr, mesh.nz
    nodes, cells = (nr+1)*nz, nr*nz
    rr, rc, rv, zr, zc, zv = [], [], [], [], [], []
    for i in range(nr+1):
        for j in range(nz):
            if i:
                rr.extend((i*nz+j, i*nz+j)); rc.extend((i*nz+j, i*nz+(j+1)%nz))
                rv.extend((1/(mesh.rf[i]*mesh.dz), -1/(mesh.rf[i]*mesh.dz)))
            if i < nr:
                zr.extend((i*nz+j, i*nz+j)); zc.extend((i*nz+j, (i+1)*nz+j))
                zv.extend((-1/mesh.radial_measure[i], 1/mesh.radial_measure[i]))
    radial = coo_matrix((rv, (rr, rc)), shape=(nodes, nodes)).tocsr()
    axial = coo_matrix((zv, (zr, zc)), shape=(cells, nodes)).tocsr()
    wr, wc, wv, vr, vc, vv = [], [], [], [], [], []
    ll, lr, rr = radial_mass(nr, mesh.dr)
    for i in range(nr):
        weight = mesh.volumes[i, 0]/MU0
        for j in range(nz):
            a, b = i*nz+j, (i+1)*nz+j
            for row, col, val in ((a, a, ll[i, 0]), (a, b, lr[i, 0]), (b, a, lr[i, 0]), (b, b, rr[i, 0])):
                wr.append(row); wc.append(col); wv.append(weight*val)
            a, b = i*nz+j, i*nz+(j+1)%nz
            for row, col, val in ((a, a, 1/3), (a, b, 1/6), (b, a, 1/6), (b, b, 1/3)):
                vr.append(row); vc.append(col); vv.append(weight*val)
    mr = coo_matrix((wv, (wr, wc)), shape=(nodes, nodes)).tocsr()
    mz = coo_matrix((vv, (vr, vc)), shape=(cells, cells)).tocsr()
    poloidal = radial.T @ mr @ radial+axial.T @ mz @ axial
    s = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]
    toroidal = diags((mesh.volumes*s/MU0).ravel())
    return block_diag((poloidal, toroidal), format='csc')


def transport_matrix(geometry, magnetic, number_density=None, resistivity=0.):
    """Return X+R: Hall skew block plus positive resistive current mass."""
    maps, weight = geometry.maps, geometry.weight.ravel()
    eta = np.broadcast_to(np.asarray(resistivity, float), geometry.r.shape).ravel()
    if not np.isfinite(eta).all() or np.any(eta < 0):
        raise ValueError('Finite nonnegative resistivity at quadrature required')
    resistive = [(p.T @ p.multiply((weight*eta)[:, None])).tocsc() for p in maps]
    if number_density is None:
        return block_diag(resistive, format='csc')
    density = np.broadcast_to(np.asarray(number_density, float), geometry.r.shape).ravel()
    if not np.isfinite(density).all() or np.any(density <= 0):
        raise ValueError('Strictly positive finite Hall carrier density required')
    field = np.asarray(magnetic).reshape(-1, 3)
    if not np.isfinite(field).all():
        raise ValueError('Finite magnetic field required')
    def block(i, j, k, sign):
        return (maps[i].T @ maps[j].multiply((sign*weight*field[:, k]/(density*E))[:, None])).tocsc()
    rt, rz, tz = block(0, 1, 2, 1), block(0, 2, 1, -1), block(1, 2, 0, 1)
    return bmat([[resistive[0], rt, rz], [-rt.T, resistive[1], tz], [-rz.T, -tz.T, resistive[2]]], format='csc')
