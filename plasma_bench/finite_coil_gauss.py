"""Unsplit Gauss evolution using finite-tube active coordinates and source fields."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0,electron_mass as ME
from scipy.sparse import block_diag,csr_matrix,eye,kron,diags
from scipy.sparse.linalg import splu
from scipy.linalg import solve,lu_factor,lu_solve
from .driven_modal_gauss import DrivenModalGauss
from .finite_coil_source import FiniteCoilMagnetic
from .regular_fluid_linear import radial_matrices
from .spitzer_wall_plasma import spitzer_samples


class FiniteCoordinates:
    def __init__(self,system,state,field):
        self.system=system;self.fluid=d=system.fluid;self.geometry=g=system.geometry;m=g.mesh
        f=d.fields(state);field.packed_increment(g)
        self.density=float(np.sum(d.volume*state[...,0])/np.sum(d.volume))
        self.pressures=np.sum(d.volume[...,None]*f['pressure'],axis=(0,1,2))/np.sum(d.volume)
        self.speed=float(np.sqrt(np.sum(d.thermal.gammas*self.pressures)/(d.mass*self.density)))
        active=np.ones(d.shape,bool);active[0,0,:,1:3]=False;active[-1,-1,:,1]=False;active[...,[0,-1],3]=False
        self.fluid_indices=np.flatnonzero(active)
        scales=self.density*np.array([1,self.speed/d.radius,self.speed/d.radius,self.speed,1,1])
        self.fluid_scale=np.broadcast_to(scales,d.shape).ravel()[self.fluid_indices]
        bscale=np.sqrt(MU0*d.mass*self.density)*self.speed
        self.psi_scale=bscale*d.radius**2;self.g_scale=bscale/d.radius
        self.npsi=m.nr*(m.nz+1);self.gshape=(g.radial_null.shape[1],g.axial_null.shape[1])
        self.nfluid=len(self.fluid_indices);self.size=self.nfluid+self.npsi+np.prod(self.gshape)
        ids=np.arange(m.nz+1,(m.nr+1)*(m.nz+1))
        psilift=csr_matrix((np.full(len(ids),self.psi_scale),(ids,np.arange(len(ids)))),shape=((m.nr+1)*(m.nz+1),len(ids)))
        self.lift=block_diag((psilift,kron(csr_matrix(g.radial_null),csr_matrix(g.axial_null))*self.g_scale),format='csr')

    def pack(self,fluid,magnetic):
        p,b=self.geometry.unpack(magnetic)
        return np.r_[np.asarray(fluid).ravel()[self.fluid_indices]/self.fluid_scale,
                     p[1:].ravel()/self.psi_scale,(self.geometry.radial_null.T@b@self.geometry.axial_null).ravel()/self.g_scale]

    def unpack(self,value):
        value=np.asarray(value,float)
        if value.shape!=(self.size,) or not np.isfinite(value).all():raise ValueError('Finite full cylinder coordinates required')
        fluid=np.zeros(self.fluid.shape);fluid.ravel()[self.fluid_indices]=value[:self.nfluid]*self.fluid_scale
        return fluid,np.asarray(self.lift@value[self.nfluid:])

    def background(self):
        d=self.fluid;state=np.zeros(d.shape);state[...,0]=self.density
        state[...,4:]=self.density*(np.log(self.pressures/(E*self.density))-(d.thermal.gammas-1)*np.log(self.density/d.thermal.n0))
        return state


class FiniteInverse:
    """Uniform fluid and resistive/inertial magnetic blocks; proposal only."""
    def __init__(self,c,state):
        self.coordinates=c;d=c.fluid;g=c.geometry;system=c.system
        f=d.fields(c.background());entropy=f['entropy'][0,0,0];h=f['specific_internal'][0,0,0]
        pressure=np.zeros(6);pressure[0]=np.sum((d.thermal.gammas-1)*h*(d.thermal.gammas-entropy));pressure[4:]=(d.thermal.gammas-1)*h
        kr=np.zeros((6,6));kr[[0,4,5],1]=[1,*entropy]
        kp=np.zeros((6,6));kp[1]=pressure/d.mass
        kz=np.zeros((6,6));kz[0,3]=1;kz[3]=pressure/d.mass;kz[4:,3]=entropy
        mass,gradient=radial_matrices(d);z=csr_matrix(-d.axial_sbp.q/d.axial_weight[:,None])
        fluid=kron(kron(mass,eye(d.nz)),kr)+kron(kron(gradient,eye(d.nz)),kp)+kron(kron(eye(mass.shape[0]),z),kz)
        fluid=fluid.tocsr()[c.fluid_indices][:,c.fluid_indices]
        self.fluid=(diags(1/c.fluid_scale)@fluid@diags(c.fluid_scale)).tocsc()
        weak=-(g.faraday.T@g.hessian@c.lift).toarray()
        current=np.vstack([solver(weak[g.offsets[k]:g.offsets[k+1]]) for k,solver in enumerate(g.solves)])
        jmetric=current.T@weak;mu=system.ion_mass*ME/(system.ion_mass+ME)
        metric=(c.lift.T@g.hessian@c.lift).toarray()+mu/(E**2*c.density)*jmetric
        if system.transport=='fixed':eta=system.resistivity
        else:eta=float(np.sum(g.weight*spitzer_samples(system.mapping,state,system.ion_mass)[2])/np.sum(g.weight))
        self.magnetic=solve(metric,-eta*jmetric,assume_a='pos');self.shift=None

    def factor(self,shift):
        if shift==self.shift:return
        self.ff=splu(eye(self.fluid.shape[0],format='csc',dtype=complex)-shift*self.fluid)
        self.mf=lu_factor(np.eye(self.magnetic.shape[0])-shift*self.magnetic);self.shift=shift

    def solve(self,value):
        n=self.coordinates.nfluid
        return np.r_[self.ff.solve(value[:n]),lu_solve(self.mf,value[n:])]

    def stages(self,value):
        k=2-np.sqrt(3);solved=self.solve(.5*value[1]-.5j/k*value[0])
        return np.stack((-2*k*solved.imag,2*solved.real))


class FiniteCoilGauss(DrivenModalGauss):
    def __init__(self,system,state,field,*,absolute_tolerance=1e-12,relative_tolerance=1e-10,max_iterations=60):
        if not np.isfinite([absolute_tolerance,relative_tolerance]).all() or min(absolute_tolerance,relative_tolerance)<=0:
            raise ValueError('Positive finite tolerances required')
        if type(max_iterations) is not int or max_iterations<1:raise ValueError('Positive iteration budget required')
        self.system=system;self.coordinates=FiniteCoordinates(system,state,field);self.inverse=FiniteInverse(self.coordinates,state)
        self.atol=absolute_tolerance;self.rtol=relative_tolerance;self.max_iterations=max_iterations
        r=np.sqrt(3)/6;self.a=np.array([[.25,.25-r],[.25+r,.25]])

    def source_value(self,field):return field.source_amplitude

    def make_field(self,packed,time):
        psi,g=self.system.geometry.unpack(packed)
        return FiniteCoilMagnetic(self.system.source,psi,g,self.system.field_value(time))
