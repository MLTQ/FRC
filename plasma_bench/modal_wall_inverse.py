"""Complex shifted plasma inverse using independent small Fourier blocks."""
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from scipy.linalg import lu_factor, lu_solve
from .modal_wall_blocks import ModalWallBlocks


class ModalWallInverse:
    def __init__(self, coordinates):
        self.coordinates = coordinates; self.blocks = ModalWallBlocks(coordinates)
        self.operators = [(self.blocks.fluid(k), *self.blocks.coupled(k)) for k in range(coordinates.nz)]
        self._shift = None

    def factor(self, shift):
        h = complex(shift)
        if not np.isfinite([h.real, h.imag]).all() or h.real <= 0: raise ValueError('Positive-real finite shift required')
        if self._shift == h: return
        factors = []
        for fluid, to_fluid, to_magnetic, magnetic in self.operators:
            f = splu(eye(fluid.shape[0], format='csc', dtype=complex)-h*fluid)
            response = f.solve(to_fluid)
            schur = np.eye(magnetic.shape[0])-h*magnetic-h*h*to_magnetic@response
            factors.append((f, response, lu_factor(schur)))
        self.factors = factors; self._shift = h

    def solve(self, value):
        if self._shift is None: raise ValueError('Factor the shifted operator first')
        c = self.coordinates; modes = c.fourier(value); result = np.empty_like(modes, complex)
        n = c.fluid_row_count; h = self._shift
        for k, ((_, _, to_magnetic, _), (factor, response, magnetic_factor)) in enumerate(zip(self.operators, self.factors)):
            fluid = factor.solve(modes[:n, k])
            magnetic = lu_solve(magnetic_factor, modes[n:, k]+h*to_magnetic@fluid)
            result[:n, k] = fluid+h*response@magnetic; result[n:, k] = magnetic
        return c.physical(result)

    def stages(self, value):
        kappa = 2-np.sqrt(3.)
        solved = self.solve(.5*value[1]-.5j/kappa*value[0])
        return np.stack((-2*kappa*solved.imag, 2*solved.real))
