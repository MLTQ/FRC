"""Consistent-mass resistive induction, actual vacuum and positive material heat."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


@dataclass
class ResistiveUpdate:
    material_a: np.ndarray
    induced: np.ndarray
    midpoint_current_over_radius: np.ndarray
    cell_heat: np.ndarray
    source_work: float
    initial_energy: float
    final_energy: float
    balance_residual: float
    linear_relative_residual: float
    vacuum_relative_residual: float


class MaterialResistivity:
    """Frozen-geometry component with scalar physical eta [ohm m], no fluid outside."""
    def __init__(self,domain,coordinates,resistivity):
        if not np.isscalar(resistivity) or not np.isfinite(resistivity) or resistivity < 0:
            raise ValueError("Finite nonnegative scalar physical resistivity required")
        self.domain,self.resistivity=domain,float(resistivity)
        self.coordinates=np.asarray(coordinates,float).copy()
        self.elements=domain.elements(self.coordinates)
        self.stiffness=domain.stiffness(self.elements)
        n=domain.fluid_node_count;nc=domain.fluid_cell_count
        batch=self.elements[0]
        self.shape=batch.field_shape[:nc]
        self.ids=batch.field_ids[:nc]
        self.weight=batch.weights[:nc]*batch.points[:nc,:,0]**3
        local=np.einsum("cq,cqi,cqj->cij",self.weight,self.shape,self.shape)
        rows=np.repeat(self.ids,4,axis=1).ravel();cols=np.tile(self.ids,(1,4)).ravel()
        self.mass_full=coo_matrix((local.ravel(),(rows,cols)),shape=self.stiffness.shape).tocsc()
        self.mass=self.mass_full[:n,:n].tocsc()
        self.mass_factor=splu(self.mass)
        self.vacuum_factor=splu(self.stiffness[n:,n:].tocsc())
        self._dt,self._factor,self._matrix=None,None,None

    def induced(self,material_a,bias=0.):
        d=self.domain;n=d.fluid_node_count
        material_a=np.asarray(material_a,float)
        if material_a.shape != (n,) or not np.all(np.isfinite(material_a)) or not np.isfinite(bias):
            raise ValueError("Finite material a coefficients and bias required")
        u=np.zeros(d.field_count);u[:n]=material_a-.5*bias
        u[n:]=self.vacuum_factor.solve(-(self.stiffness[n:,:n]@u[:n]))
        return u

    def current_over_radius(self,material_a,bias=0.):
        u=self.induced(material_a,bias)
        return self.mass_factor.solve((self.stiffness@u)[:self.domain.fluid_node_count])/MU0

    def current_quadrature(self,coefficients):
        r=self.elements[0].points[:self.domain.fluid_cell_count,:,0]
        return r*np.einsum("ci,cqi->cq",coefficients[self.ids],self.shape)

    def step(self,material_a,dt,*,bias_start=0.,bias_end=None):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Finite positive resistive step required")
        bias_end=bias_start if bias_end is None else bias_end
        if not np.isfinite(bias_start) or not np.isfinite(bias_end):
            raise ValueError("Finite prescribed bias endpoints required")
        n=self.domain.fluid_node_count;k=self.stiffness;m=self.mass_full
        u0=self.induced(material_a,bias_start);delta=float(bias_end-bias_start)
        source=np.zeros(self.domain.field_count);source[:n]=.5*delta
        if self.resistivity == 0:
            u1=self.induced(material_a,bias_end);linear_residual=0.
        else:
            alpha=.5*dt*self.resistivity/MU0
            if self._dt != dt:
                self._matrix=(m+alpha*k).tocsc()
                self._factor=splu(self._matrix);self._dt=dt
            rhs=m@(u0-source)-alpha*(k@u0)
            rhs[n:]=0.  # Exact stationary vacuum row, not a fluid update.
            u1=self._factor.solve(rhs)
            linear_residual=float(np.linalg.norm(self._matrix@u1-rhs)/max(np.linalg.norm(rhs),np.finfo(float).tiny))
        midpoint=.5*(u0+u1);ku=k@midpoint
        q=self.mass_factor.solve(ku[:n])/MU0
        qfield=np.einsum("ci,cqi->cq",q[self.ids],self.shape)
        heat=2*np.pi*self.resistivity*dt*np.sum(self.weight*qfield*qfield,axis=1)
        work=float(-np.pi/MU0*delta*np.sum(ku[:n]))
        energy0=float(np.pi/MU0*np.dot(u0,k@u0));energy1=float(np.pi/MU0*np.dot(u1,k@u1))
        load=k@u1
        vacuum=float(np.linalg.norm(load[n:])/max(np.linalg.norm(load[:n]),np.finfo(float).tiny))
        if not np.all(np.isfinite(u1)) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(heat)) or not np.isfinite(energy0+energy1+work):
            raise ValueError("Nonfinite resistive proposal; no state is committed")
        if linear_residual > 1e-9 or vacuum > 1e-8:
            raise ValueError("Unconverged material/vacuum resistive solve")
        material_new=np.asarray(material_a,float).copy() if self.resistivity == 0 else u1[:n]+.5*bias_end
        return ResistiveUpdate(material_new,u1,q,heat,work,energy0,energy1,
                               float(energy1-energy0+heat.sum()-work),linear_residual,vacuum)
