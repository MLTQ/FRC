"""Unchanged full Gauss residual solve with a scalable uniform modal inverse."""
import numpy as np
from .implicit_plasma_gauss import ImplicitPlasmaGauss
from .modal_wall_coordinates import ModalWallCoordinates
from .modal_wall_inverse import ModalWallInverse


class ModalWallGauss(ImplicitPlasmaGauss):
    def __init__(self, system, state, field, *, absolute_tolerance=1e-12,
                 relative_tolerance=1e-10, max_iterations=60):
        if not np.isfinite([absolute_tolerance, relative_tolerance]).all() or min(absolute_tolerance, relative_tolerance) <= 0:
            raise ValueError('Positive finite nonlinear residual tolerances required')
        if type(max_iterations) is not int or max_iterations < 1: raise ValueError('Positive iteration budget required')
        self.system = system; self.coordinates = ModalWallCoordinates(system, state, field)
        self.inverse = ModalWallInverse(self.coordinates)
        self.atol = absolute_tolerance; self.rtol = relative_tolerance; self.max_iterations = max_iterations
        r = np.sqrt(3.)/6; self.a = np.array([[.25, .25-r], [.25+r, .25]])
