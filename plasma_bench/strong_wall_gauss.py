"""Full nonlinear Gauss equations in a boundary-conforming coordinate space."""
import numpy as np
from .implicit_plasma_gauss import ImplicitPlasmaGauss
from .strong_wall_coordinates import StrongWallCoordinates
from .strong_wall_preconditioner import StrongWallPreconditioner
from .complex_plasma_inverse import ComplexPlasmaInverse


class StrongWallGauss(ImplicitPlasmaGauss):
    def __init__(self, system, state, field, *, absolute_tolerance=1e-12,
                 relative_tolerance=1e-10, max_iterations=60):
        if not np.isfinite([absolute_tolerance, relative_tolerance]).all() or min(absolute_tolerance, relative_tolerance) <= 0:
            raise ValueError('Positive finite nonlinear residual tolerances required')
        if type(max_iterations) is not int or max_iterations < 1:
            raise ValueError('Positive nonlinear iteration budget required')
        self.system = system; self.coordinates = StrongWallCoordinates(system, state, field)
        self.preconditioner = StrongWallPreconditioner(self.coordinates)
        self.inverse = ComplexPlasmaInverse(self.preconditioner)
        self.atol = absolute_tolerance; self.rtol = relative_tolerance; self.max_iterations = max_iterations
        r = np.sqrt(3.)/6
        self.a = np.array([[.25, .25-r], [.25+r, .25]])
