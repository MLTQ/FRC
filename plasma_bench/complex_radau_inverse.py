"""Conjugate shifted inverse for real two-stage Radau IIA."""
import numpy as np
from .complex_plasma_inverse import ComplexPlasmaInverse


class ComplexRadauInverse(ComplexPlasmaInverse):
    def stages(self, value):
        a = 1/9; b = 2*np.sqrt(2.)/9
        solved = self.solve(.5*value[1]+.5j*(a*value[1]-value[0])/b)
        return np.stack((2*(a*solved.real-b*solved.imag), 2*solved.real))
