"""Full stationary coupled Jacobian in the physical perturbation energy metric."""
import time
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.linalg import block_diag, cholesky, solve_triangular
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_plasma_step import field_from_packed
from .regular_plasma_reference import initialize


def run_case(spec):
    start = time.monotonic()
    ref, d, u, system, field = initialize(dict(spec, kind='wave'))
    n = ref.number_density; cs = ref.sound_speed
    u[...] = 0.; u[..., 0] = n
    pressures = np.array([ref.electron_pressure, ref.ion_pressure])
    entropy = np.log(pressures/(E*n))-(d.thermal.gammas-1)*np.log(n/d.thermal.n0)
    u[..., 4:] = n*entropy
    zero = np.zeros(system.geometry.hessian.shape[0])
    field = field_from_packed(system.geometry, zero, ref.bias)
    active = np.ones(d.shape, bool); active[0, 0, :, 1:3] = False
    magnetic_mask = np.ones(zero.shape, bool)
    nz = system.geometry.mesh.nz; nr = system.geometry.mesh.nr
    magnetic_mask[:nz] = False; magnetic_mask[nr*nz:(nr+1)*nz] = False
    scales = np.broadcast_to(np.array([n, n*cs/d.radius, n*cs/d.radius, n*cs, n, n]), d.shape)[active]
    bscale = np.sqrt(MU0*n*d.mass)*cs
    magnetic_scales = np.r_[np.full((nr+1)*nz, bscale*d.radius**2), np.full(nr*nz, bscale/d.radius)][magnetic_mask]
    scale = np.r_[scales, magnetic_scales]
    nf = len(scales); size = len(scale); epsilon = spec.get('epsilon', 1e-5)

    def action(y):
        state = u.copy(); state[active] += y[:nf]*scales
        b = zero.copy(); b[magnetic_mask] = y[nf:]*magnetic_scales
        rate = system.derivative(state, field_from_packed(system.geometry, b, ref.bias))
        return np.r_[rate.fluid[active], rate.magnetic[magnetic_mask]]/scale

    columns = []
    for k in range(size):
        y = np.zeros(size); y[k] = epsilon
        columns.append((action(y)-action(-y))/(2*epsilon))
    matrix = np.column_stack(columns)
    # Independently differentiated physical total energy at uniform rest.
    f = d.fields(u); h = f['specific_internal']; gammas = d.thermal.gammas
    blocks = np.zeros(d.shape[:-1]+(6, 6))
    blocks[..., 0, 0] = np.sum(h*((gammas-entropy)*(gammas-1-entropy)+entropy), axis=-1)/n
    for k in range(2):
        blocks[..., 0, 4+k] = blocks[..., 4+k, 0] = h[..., k]*(gammas[k]-1-entropy[k])/n
        blocks[..., 4+k, 4+k] = h[..., k]/n
    blocks[..., 1, 1] = blocks[..., 2, 2] = d.mass*d.r**2/n
    blocks[..., 3, 3] = d.mass/n
    blocks *= d.volume[..., None, None]
    fluid_metric = block_diag(*blocks.reshape(-1, 6, 6))[np.ix_(active.ravel(), active.ravel())]
    model = RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)
    indices = np.flatnonzero(magnetic_mask)
    ampere = np.column_stack([model.current(np.eye(len(zero))[i]) for i in indices])
    inertia_ampere = np.column_stack([model.inertia @ a for a in ampere.T])
    magnetic_metric = model.hessian[np.ix_(indices, indices)].toarray()+ampere.T@inertia_ampere
    metric = block_diag(fluid_metric, magnetic_metric)*scale[:, None]*scale[None, :]
    root = cholesky(metric, lower=False)
    weighted = solve_triangular(root.T, (root@matrix).T, lower=True).T
    eigenvalues = np.linalg.eigvals(weighted)
    radius = float(np.max(abs(eigenvalues)))
    result = dict(spec=spec, status='complete', size=size, stationary_residual=float(np.linalg.norm(action(np.zeros(size)))),
        spectral_radius=radius, maximum_real=float(eigenvalues.real.max()),
        maximum_real_over_radius=float(eigenvalues.real.max()/radius),
        energy_symmetric_defect=float(np.linalg.norm(weighted+weighted.T)/np.linalg.norm(weighted)),
        imaginary_axis_rk4_limit=2*np.sqrt(2)/radius,
        electron_cyclotron_frequency=abs(E*ref.bias/ME), elapsed_seconds=time.monotonic()-start)
    return result, dict(operator=matrix, energy_metric=metric, energy_weighted_operator=weighted,
                        eigenvalues=eigenvalues, scale=scale)
