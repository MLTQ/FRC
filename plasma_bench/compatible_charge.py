"""Divergence-free dual charge flux paired with the high-order magnetic Hessian."""
import numpy as np
from .dual_current import DualCurrent


class CompatibleCharge:
    """The graph is unchanged; its stream uses the full physical magnetic metric."""
    def __init__(self, geometry):
        dual=DualCurrent(geometry.mesh)
        self.mesh=dual.mesh;self.dx=dual.dx;self.volumes=dual.volumes
        self.r=dual.r;self.z=dual.z
        self.scalar=dual.scalar;self.gradient=dual.gradient
        self.divergence=dual.divergence;self.faraday=dual.faraday
        self.field_geometry=geometry
        nodes=(self.mesh.nr+1)*self.mesh.nz
        self.toroidal_hessian=geometry.hessian[nodes:,nodes:]

    def current_flux(self, g, wall_btheta=0.):
        g=np.asarray(g,float)
        if g.shape!=(self.mesh.nr,self.mesh.nz) or not np.isfinite(g).all():
            raise ValueError('Finite toroidal cell coefficients required')
        if not np.isfinite(wall_btheta) or wall_btheta!=0:
            raise ValueError('Compatible charge currently requires zero natural wall Btheta')
        h=(self.toroidal_hessian@g.ravel()).reshape(g.shape)
        stream=h/(self.dx[:,None]*self.mesh.dz)
        radial=np.zeros((self.mesh.nr+2,self.mesh.nz))
        radial[1:-1]=np.roll(stream,1,axis=1)-stream
        axial=np.diff(np.concatenate((np.zeros((1,self.mesh.nz)),stream,
                                     np.zeros((1,self.mesh.nz))),axis=0),axis=0)
        return radial,axial

    def work(self, g, wall_btheta, radial_coefficient, axial_electric):
        radial,axial=self.current_flux(g,wall_btheta)
        work=np.sum(radial[1:-1]*self.dx[:,None]*radial_coefficient)
        work+=np.sum(axial*self.mesh.dz*axial_electric)
        rate=self.faraday(radial_coefficient,axial_electric)
        magnetic=g.ravel()@(self.toroidal_hessian@rate.ravel())
        return dict(current_work_w=float(work),boundary_magnetic_input_w=0.,
                    magnetic_rate_w=float(magnetic))
