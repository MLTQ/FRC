"""Fixed-material linear current response with physical electron inertia and vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from scipy.sparse import bmat, coo_matrix, diags
from scipy.sparse.linalg import splu
from .material_resistivity import MaterialResistivity


@dataclass
class CurrentUpdate:
    material_a: np.ndarray
    induced: np.ndarray
    current_over_radius: np.ndarray
    midpoint_current_over_radius: np.ndarray
    cell_heat: np.ndarray
    source_work: float
    initial_magnetic_energy: float
    final_magnetic_energy: float
    initial_current_energy: float
    final_current_energy: float
    balance_residual: float
    ampere_relative_residual: float
    ohm_relative_residual: float
    vacuum_relative_residual: float


class MaterialCurrent:
    """Etheta=eta Jtheta + me/(ne e²) d_t Jtheta on fixed geometry.

    Coefficients may be scalars, cell arrays or callbacks at physical quadrature
    points. Density is positive inside; exterior has no carriers. Density=None
    explicitly selects the massless comparison, not a physical density floor.
    """
    def __init__(self,domain,coordinates,*,electron_density,resistivity):
        self.base=MaterialResistivity(domain,coordinates,0.)
        self.domain=domain;self.elements=self.base.elements
        self.stiffness=self.base.stiffness;self.mass=self.base.mass
        self.weight=self.base.weight;self.shape=self.base.shape;self.ids=self.base.ids
        self.resistivity=self._sample(resistivity,"resistivity",positive=False)
        self.electron_density=None if electron_density is None else self._sample(electron_density,"electron density",positive=True)
        self.de_squared=np.zeros_like(self.weight) if electron_density is None else ME/(MU0*E**2*self.electron_density)
        self.inertia=self._mass(self.de_squared)
        self.diffusion=self._mass(self.resistivity/MU0)
        self._dt=None;self._factor=None

    def _sample(self,value,name,*,positive):
        points=self.elements[0].points[:self.domain.fluid_cell_count]
        raw=np.asarray(value(points[...,0],points[...,1]) if callable(value) else value,float)
        if raw.shape==(self.domain.fluid_cell_count,):raw=raw[:,None]
        try:result=np.broadcast_to(raw,self.weight.shape).copy()
        except ValueError as error:raise ValueError(f"Invalid {name} shape") from error
        if not np.all(np.isfinite(result)) or np.any(result<=0 if positive else result<0):
            raise ValueError(f"Finite {'positive' if positive else 'nonnegative'} {name} required inside plasma")
        return result

    def _mass(self,coefficient):
        local=np.einsum("cq,cqi,cqj->cij",self.weight*coefficient,self.shape,self.shape)
        rows=np.repeat(self.ids,4,axis=1).ravel();cols=np.tile(self.ids,(1,4)).ravel()
        n=self.domain.fluid_node_count
        return coo_matrix((local.ravel(),(rows,cols)),shape=(n,n)).tocsc()

    def induced(self,material_a,bias=0.):
        return self.base.induced(material_a,bias)

    def current_over_radius(self,material_a,bias=0.):
        return self.base.current_over_radius(material_a,bias)

    def current_quadrature(self,q):
        return self.base.current_quadrature(q)

    def energies(self,u,w):
        return (float(np.pi/MU0*np.dot(u,self.stiffness@u)),
                float(np.pi/MU0*np.dot(w,self.inertia@w)))

    def step(self,material_a,dt,*,bias_start=0.,bias_end=None):
        if not np.isfinite(dt) or dt<=0:raise ValueError("Finite positive current step required")
        bias_end=bias_start if bias_end is None else bias_end
        if not np.isfinite(bias_start) or not np.isfinite(bias_end):
            raise ValueError("Finite prescribed bias endpoints required")
        if self.electron_density is None and not np.any(self.resistivity):
            ideal=self.base.step(material_a,dt,bias_start=bias_start,bias_end=bias_end)
            q=self.base.current_over_radius(ideal.material_a,bias_end)
            load=self.stiffness@ideal.induced
            current_load=self.base.mass_full[:,:self.domain.fluid_node_count]@(MU0*q)
            ampere=float(np.linalg.norm(load-current_load)/max(np.linalg.norm(load),np.linalg.norm(current_load),np.finfo(float).tiny))
            return CurrentUpdate(ideal.material_a,ideal.induced,q,ideal.midpoint_current_over_radius,
                                 ideal.cell_heat,ideal.source_work,ideal.initial_energy,ideal.final_energy,
                                 0.,0.,ideal.balance_residual,ampere,0.,ideal.vacuum_relative_residual)
        u0=self.induced(material_a,bias_start)
        n=self.domain.fluid_node_count;nf=self.domain.field_count
        w0=self.base.mass_factor.solve((self.stiffness@u0)[:n])
        m=self.base.mass_full[:n,:];h=self.inertia;r=self.diffusion
        if self._dt!=dt:
            matrix=bmat([[self.stiffness,-self.base.mass_full[:,:n]],
                         [m,h+.5*dt*r]],format="csc")
            current_diagonal=(h+.5*dt*r).diagonal()
            current_diagonal=np.where(current_diagonal>0,current_diagonal,self.mass.diagonal()*self.domain.radius**2)
            self._scale=1/np.sqrt(np.r_[self.stiffness.diagonal(),current_diagonal])
            scaling=diags(self._scale,format="csc")
            self._factor=splu(scaling@matrix@scaling);self._dt=dt
        source=np.zeros(nf);source[:n]=.5*(bias_end-bias_start)
        rhs_ohm=-(m@source)-dt*(r@w0)
        proposal=self._scale*self._factor.solve(self._scale*np.r_[np.zeros(nf),rhs_ohm])
        du,dw=proposal[:nf],proposal[nf:]
        u1,w1=u0+du,w0+dw;wm=w0+.5*dw
        load=self.stiffness@u1;current_load=self.base.mass_full[:,:n]@w1
        ampere=float(np.linalg.norm(load-current_load)/max(np.linalg.norm(load),np.linalg.norm(current_load),np.finfo(float).tiny))
        ohm_terms=(m@du,h@dw,.5*dt*(r@dw))
        ohm_scale=sum(np.linalg.norm(t) for t in ohm_terms)+np.linalg.norm(rhs_ohm)
        ohm=float(np.linalg.norm(sum(ohm_terms)-rhs_ohm)/max(ohm_scale,np.finfo(float).tiny))
        vacuum=float(np.linalg.norm(load[n:])/max(np.linalg.norm(load[:n]),np.finfo(float).tiny))
        samples=np.einsum("ci,cqi->cq",wm[self.ids],self.shape)
        heat=2*np.pi*dt/MU0**2*np.sum(self.weight*self.resistivity*samples**2,axis=1)
        work=float(-2*np.pi/MU0*np.dot(wm,m@source))
        b0,e0=self.energies(u0,w0);b1,e1=self.energies(u1,w1)
        if not np.all(np.isfinite(proposal)) or not np.all(np.isfinite(heat)) or not np.isfinite(b0+e0+b1+e1+work):
            raise ValueError("Nonfinite current proposal; no accepted state changed")
        if ampere>1e-8 or ohm>1e-8 or vacuum>1e-8:
            raise ValueError(f"Unconverged current/vacuum solve: {ampere:g}, {ohm:g}, {vacuum:g}")
        material_new=np.asarray(material_a,float)+du[:n]+.5*(bias_end-bias_start)
        return CurrentUpdate(material_new,u1,w1/MU0,wm/MU0,heat,work,
                             b0,b1,e0,e1,float(b1-b0+e1-e0+heat.sum()-work),ampere,ohm,vacuum)
