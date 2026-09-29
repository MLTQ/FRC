"""Explicit smooth-transport comparison of the complete finite-mass model."""
from .finite_mass_mhd import FiniteMassMHD
from .smooth_finite_mass_rhs import rhs


class SmoothFiniteMassMHD(FiniteMassMHD):
    def __init__(self,*args,transport_quadrature=5,entropy_transport='smooth',**kwargs):
        if type(transport_quadrature) is not int or transport_quadrature<5:
            raise ValueError('At least fifth-order transport quadrature required')
        if entropy_transport not in ('smooth','limited'):raise ValueError('Explicit smooth or limited entropy transport required')
        self.transport_quadrature=transport_quadrature;self.entropy_transport=entropy_transport
        super().__init__(*args,**kwargs)

    def derivative(self,state,time):
        return rhs(self.mesh,state,self.parameters,time=time,source=self.source,quadrature=self.transport_quadrature,entropy_transport=self.entropy_transport)
