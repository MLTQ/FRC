"""Compatible discrete vacuum fields for known external coils, before adding any plasma response."""
import numpy as np
from .hybrid3d import curl_edges, curl_faces


class VacuumCoilBasis:
    """Make the known coil field current-free in gas, retaining its boundary-normal flux.

    This is a source-only magnetostatic solve. It never projects the evolving
    plasma field or restricts an azimuthal mode. Poloidal coil sources have no
    circulation around the annular gas region's axial hole.
    """
    def __init__(self, box, gas):
        from .hybrid3d_device import GapProjector, GAP, CONDUCTOR
        self.box = box
        self.project = GapProjector(np.where(gas, GAP, CONDUCTOR), box.d)
        self.laplacian = sum((4*np.sin(np.pi*np.fft.fftfreq(box.shape[a]))**2/box.d[a]**2)
                             .reshape(tuple(box.shape[a] if q == a else 1 for q in range(3))) for a in range(3))
        self.laplacian[0, 0, 0] = 1.

    def __call__(self, edges):
        original = np.asarray(curl_edges(*edges, self.box.d))
        target = np.asarray(self.project([b.copy() for b in original]))
        delta = target-original
        # delta is solenoidal and compact inside gas. The discrete Coulomb
        # inverse restores an edge A with curl A = projected source B, so
        # source induction and reciprocity use the very same field basis.
        correction = np.array([np.fft.ifftn(np.fft.fftn(j)/self.laplacian).real
                               for j in curl_faces(*delta, self.box.d)])
        return np.asarray(edges)+correction
