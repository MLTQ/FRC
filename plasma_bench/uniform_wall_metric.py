"""SPD uniform-density current-inertia preconditioner for the actual wall solve."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.linalg import block_diag, cho_factor, cho_solve, null_space
from .compatible_geometry import gram


class UniformWallMetric:
    def __init__(self, model, normalized_lift):
        geo = model.geometry; mesh = model.mesh; nr, nz = mesh.nr, mesh.nz
        self.nz = nz; radial = null_space(geo.wall_trace[None, :])
        lift = block_diag(np.eye(nr-1), radial)
        # The old Krylov lift has orthonormal unscaled columns. Its Euclidean
        # column norms recover positive physical Hessian normalization scales.
        actual_scale = np.sqrt(np.asarray(normalized_lift.power(2).sum(axis=0)).ravel()).reshape(-1, nz)
        row_scale = actual_scale.mean(axis=1)
        self.ratio = row_scale[:, None]/actual_scale
        scaled_lift = lift*row_scale[None, :]
        density = float(np.sum(geo.weight*model.density)/np.sum(geo.weight))
        self.reference_density = density
        coefficient = model.reduced_mass/(E**2*density)
        zn = np.fft.fft(gram(geo.Nz, geo.wz)[:, 0].toarray().ravel()).real
        ze = np.fft.fft(gram(geo.Ez, geo.wz)[:, 0].toarray().ravel()).real
        zd = np.fft.fft(gram(geo.dNz, geo.wz)[:, 0].toarray().ravel()).real
        mr = gram(geo.radial_e, geo.wr).toarray(); mt = gram(geo.theta_e, geo.wr).toarray()
        mz = gram(geo.Nx, geo.wr).toarray()
        hr = gram(geo.radial_b[:, 1:-1], geo.wr).toarray()
        hz = gram(geo.dNx[:, 1:-1], geo.wr).toarray()
        self.factors = []
        for k in range(nz):
            mass = block_diag(mr*zn[k], mt*zn[k], mz*ze[k])
            hessian = block_diag((hr*zd[k]+hz*zn[k])/MU0, mr*ze[k]/MU0)
            curl = np.zeros((2*nr-1, 3*nr), complex)
            curl[:nr-1, nr:2*nr-1] = -np.diag(mesh.rf[1:-1])
            curl[nr-1:, :nr] = np.eye(nr)*(1-np.exp(2j*np.pi*k/nz))/mesh.dz
            curl[nr-1:, 2*nr-1:] = np.diff(np.eye(nr+1), axis=0)/mesh.radial_measure[:, None]
            current = -np.linalg.solve(mass, curl.conj().T@hessian@scaled_lift)
            metric = scaled_lift.conj().T@hessian@scaled_lift+coefficient*current.conj().T@mass@current
            # Only the approximate inverse is symmetrized at floating roundoff.
            # The accepted physical operator/residual uses original matrices.
            self.factors.append(cho_factor(.5*(metric+metric.conj().T), lower=True))

    def solve(self, value):
        value = np.asarray(value, float)
        if value.size != self.ratio.size or not np.isfinite(value).all():
            raise ValueError('Finite complete normalized Krylov vector required')
        modes = np.fft.fft(self.ratio*value.reshape(self.ratio.shape), axis=1)
        result = np.empty_like(modes)
        for k, factor in enumerate(self.factors): result[:, k] = cho_solve(factor, modes[:, k])
        return (self.ratio*np.fft.ifft(result, axis=1).real).ravel()
