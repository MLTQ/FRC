"""Read-only n=0 flux and azimuthal density diagnostics for saved 3D hybrid states."""
import numpy as np
from examples.nucte_compare import closed_region


def poloidal_average(box, device, bz):
    """Return (radial nodes, axial nodes, psi, region) from annular face averages.

    psi(r,z) = integral_0^r <Bz> r dr. This is the axisymmetric component,
    not a separatrix of the full nonaxisymmetric field. Bz lives on z faces.
    """
    lo, d = np.asarray(box.lo), box.d
    xh = [lo[a]+d[a]*(np.arange(box.shape[a])+.5) for a in (0, 1)]
    x, y = np.meshgrid(*xh, indexing="ij")
    edges = np.arange(0, device.coil_radius+d[0], d[0])
    bins = np.digitize(np.hypot(x, y), edges)-1
    mean = np.array([bz[bins == q].mean(axis=0) if (bins == q).any()
                     else np.zeros(box.shape[2]) for q in range(len(edges)-1)])
    rc = .5*(edges[:-1]+edges[1:])
    psi = np.concatenate((np.zeros((1, box.shape[2])),
                          np.cumsum(mean*rc[:, None]*d[0], axis=0)))
    z = lo[2]+d[2]*np.arange(box.shape[2])
    gas_z = np.abs(z) < device.tube_half_length
    wall = int(np.searchsorted(edges, device.wall_radius))
    r, z, psi = edges[:wall+1], z[gas_z], psi[:wall+1, gas_z]
    return r, z, psi, closed_region(r, z, psi)


def density_modes(box, device, density):
    """Midplane |sum(n exp(-im theta))|/sum(n), in 0.3..0.9 wall radii.

    These are instantaneous amplitudes, including particle noise; they are
    not fitted instability growth rates.
    """
    lo, d = np.asarray(box.lo), box.d
    xs = [lo[a]+d[a]*np.arange(box.shape[a]) for a in (0, 1)]
    x, y = np.meshgrid(*xs, indexing="ij")
    r, phi = np.hypot(x, y), np.arctan2(y, x)
    ring = (r > .3*device.wall_radius) & (r < .9*device.wall_radius)
    n = density[:, :, box.shape[2]//2][ring]
    m0 = max(float(abs(n.sum())), 1e-300)
    return {f"n{m}": float(abs(np.sum(n*np.exp(-1j*m*phi[ring])))/m0)
            for m in (1, 2, 3)}
