"""Complex shifted block inverse for real Gauss collocation stages."""
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import splu
from scipy.linalg import lu_factor, lu_solve


class ComplexPlasmaInverse:
    def __init__(self, preconditioner):
        self.operator = preconditioner; self._shift = None

    def factor(self, shift):
        h = complex(shift)
        if not np.isfinite([h.real, h.imag]).all() or h.real <= 0:
            raise ValueError('Finite shift with positive real part required')
        if self._shift != h:
            p = self.operator
            self.fluid_factor = splu(eye(p.fluid.shape[0], format='csc', dtype=complex)-h*p.fluid)
            self.fluid_response = self.fluid_factor.solve(p.to_fluid.astype(complex))
            schur = np.eye(p.magnetic.shape[0], dtype=complex)-h*p.magnetic-h*h*p.to_magnetic@self.fluid_response
            self.magnetic_factor = lu_factor(schur)
            self._shift = h

    def solve(self, value):
        if self._shift is None: raise ValueError('Factor the shifted operator first')
        p = self.operator; n = p.coordinates.nfluid; h = self._shift
        fluid = self.fluid_factor.solve(np.asarray(value[:n], complex))
        magnetic = lu_solve(self.magnetic_factor, value[n:]+h*p.to_magnetic@fluid)
        return np.r_[fluid+h*self.fluid_response@magnetic, magnetic]

    def stages(self, value):
        """Apply (I-dt A_Gauss tensor L)^-1 to two real stage residuals."""
        kappa = 2-np.sqrt(3.)
        solved = self.solve(.5*value[1]-.5j/kappa*value[0])
        return np.stack((-2*kappa*solved.imag, 2*solved.real))
