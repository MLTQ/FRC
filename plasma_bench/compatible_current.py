"""Constant-density physical-mass current/inertia on compatible high-order fields."""
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME
from scipy.sparse import block_diag
from .compatible_geometry import CompatibleGeometry


class CompatibleCurrent:
    def __init__(self,mesh,number_density,ion_mass,*,degree=5,quadrature_order=10):
        if not np.isfinite(number_density) or number_density<=0 or not np.isfinite(ion_mass) or ion_mass<=0:
            raise ValueError('Positive finite physical carrier density and ion mass required')
        self.mesh=mesh;self.geometry=CompatibleGeometry(mesh,degree=degree,quadrature_order=quadrature_order)
        self.number_density=number_density;self.ion_mass=ion_mass;self.electron_mass=ME
        self.alpha=ME/(ion_mass+ME);self.beta=ion_mass/(ion_mass+ME)
        self.delta=self.beta-self.alpha;self.reduced_mass=ion_mass*self.alpha
        self.mass=block_diag(self.geometry.mass,format='csc');self.hessian=self.geometry.hessian
        self.inertia=self.reduced_mass/(number_density*E*E)*self.mass

    def current(self,magnetic):
        weak=-self.geometry.faraday.T@(self.hessian@magnetic)
        return np.concatenate([solve(v) for solve,v in zip(self.geometry.solves,self.geometry.split_electric(weak))])
