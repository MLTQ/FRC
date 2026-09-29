"""Massless infinite magnetic vacuum and analytic geometric virtual work."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse.linalg import splu


@dataclass
class MovingMagneticState:
    elements: tuple
    stiffness: object
    induced: np.ndarray
    load: np.ndarray
    energy: float
    geometry_gradient: np.ndarray
    vacuum_residual: float
    physical_energy: float


class MovingField:
    def __init__(self, domain):
        self.domain = domain

    def solve(self, coordinates, active_induced):
        d = self.domain; n = d.fluid_node_count
        active_induced = np.asarray(active_induced,float)
        if active_induced.shape != (n,) or not np.all(np.isfinite(active_induced)):
            raise ValueError("Finite induced coefficients required at material nodes")
        elements = d.elements(coordinates)
        matrix = d.stiffness(elements)
        u = np.zeros(d.field_count);u[:n] = active_induced
        u[n:] = splu(matrix[n:,n:].tocsc()).solve(-(matrix[n:,:n]@active_induced))
        load = matrix@u
        gradient = np.zeros((d.node_count,2))
        physical = 0.
        for batch in elements:
            local = batch.shape_gradient(u)
            np.add.at(gradient,batch.geometry_ids.ravel(),local.reshape(-1,2))
            br,bz = self.fields(batch,u)
            physical += np.sum(batch.weights*batch.points[...,0]*(br*br+bz*bz))
        physical += np.dot(d.physical_exterior_diagonal,u[d.node_count:]**2)
        scale = np.pi/MU0
        residual = np.linalg.norm(load[n:])/max(np.linalg.norm(load[:n]),np.finfo(float).tiny)
        return MovingMagneticState(elements,matrix,u,2*scale*load,
                                   float(scale*np.dot(u,load)), scale*(d.motion.T@gradient),
                                   float(residual),float(scale*physical))

    @staticmethod
    def fields(elements, coefficients, bias=0.):
        a = np.einsum("cn,cqn->cq",coefficients[elements.field_ids],elements.field_shape)
        gradient = np.einsum("cn,cqna->cqa",coefficients[elements.field_ids],elements.field_gradient)
        r = elements.points[...,0]
        return -r*gradient[...,1],2*a+r*gradient[...,0]+bias
