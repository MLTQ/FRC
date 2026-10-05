"""Third-order L-stable Radau IIA for the actual coupled plasma equations."""
import numpy as np
from .implicit_plasma_gauss import ImplicitPlasmaGauss
from .complex_radau_inverse import ComplexRadauInverse
from .implicit_collocation_step import step


class ImplicitPlasmaRadau(ImplicitPlasmaGauss):
    def __init__(self, system, state, field, **options):
        super().__init__(system, state, field, **options)
        self.a = np.array([[5/12, -1/12], [3/4, 1/4]])
        self.weights = self.a[1].copy(); self.endpoint_stage = 1
        self.shift = 1/3+1j*np.sqrt(2.)/6
        self.inverse = ComplexRadauInverse(self.preconditioner)

    def step(self, state, field, dt):
        return step(self, state, field, dt)
