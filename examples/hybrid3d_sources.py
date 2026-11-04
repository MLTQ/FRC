"""Axisymmetric coil flux mapped once to 3D Yee edges; saved bank currents prescribe its time dependence."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.constants import mu_0 as MU0
from plasma_bench.hybrid3d import curl_edges, curl_faces
from plasma_bench.hybrid3d_coils import VacuumCoilBasis


def flux_edges(box, grid, psi):
    """A=(-y,x,0) psi/r² on Yee edges; its discrete curl is solenoidal.

    Only the initial/coil fields use axisymmetric flux. Evolving fields and
    plasma responses are never averaged or reconstructed from this function.
    """
    scale = float(np.max(np.abs(psi))) or 1.
    interp = RegularGridInterpolator((grid.r, grid.z), psi/scale, method='cubic',
                                     bounds_error=False, fill_value=None,
                                     solver_args=dict(atol=1e-13, rtol=1e-10))
    axis = RegularGridInterpolator((grid.z,), (16*psi[1]-psi[2])/(12*grid.dr**2),
                                   bounds_error=False, fill_value=None)
    xs = [box.lo[a]+box.d[a]*np.arange(box.shape[a]) for a in range(3)]
    xh = [x+.5*box.d[a] for a, x in enumerate(xs)]

    def factor(x, y, z):
        rt = np.hypot(x, y)
        r, zz = np.minimum(rt, grid.r[-1]), np.clip(z, grid.z[0], grid.z[-1])
        val = scale*interp(np.column_stack((r.ravel(), zz.ravel()))).reshape(r.shape)
        out = val/np.maximum(rt, 1e-30)**2
        small = r < .5*grid.dr
        if small.any():
            out[small] = axis(zz[small, None]).ravel()
        return out
    x, y, z = np.meshgrid(xh[0], xs[1], xs[2], indexing='ij')
    ax = -y*factor(x, y, z)
    x, y, z = np.meshgrid(xs[0], xh[1], xs[2], indexing='ij')
    ay = x*factor(x, y, z)
    return [ax, ay, np.zeros(box.shape)]


class BankWaveform:
    """Prescribed vacuum coil-current history; no replay of the 2D plasma response."""
    def __init__(self, box, grid, unit_flux, history, start_s, gas_cells=None):
        project = VacuumCoilBasis(box, gas_cells) if gas_cells is not None else lambda a: a
        self.edges = np.array([project(flux_edges(box, grid, psi)) for psi in unit_flux])
        self.faces = np.array([curl_edges(*a, box.d) for a in self.edges])
        self.current_basis = np.array([np.asarray(curl_faces(*b, box.d))/MU0 for b in self.faces])
        self.times = np.array([h['time_s']-start_s for h in history])
        self.currents = np.array([h['circuit']['vacuum_current_a'] for h in history])
        if np.any(np.diff(self.times) <= 0):
            raise ValueError('Bank waveform times must increase strictly')

    def coefficients(self, time_s):
        if time_s < self.times[0]-1e-12 or time_s > self.times[-1]+1e-12:
            raise ValueError('3D continuation exceeds saved coil-current waveform')
        return np.array([np.interp(time_s, self.times, self.currents[:, a])
                         for a in range(self.currents.shape[1])])

    def __call__(self, time_s):
        return np.einsum('s,saxyz->axyz', self.coefficients(time_s), self.faces)

    def edge_derivative(self, time_s):
        k = int(np.clip(np.searchsorted(self.times, time_s, side='right')-1, 0, len(self.times)-2))
        slope = (self.currents[k+1]-self.currents[k])/(self.times[k+1]-self.times[k])
        return np.einsum('s,saxyz->axyz', slope, self.edges)

    def edge_current(self, time_s):
        return np.einsum('s,saxyz->axyz', self.coefficients(time_s), self.current_basis)
