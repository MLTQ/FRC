"""Full-vector material ideal magnetic stress, torque and current reconstruction."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse.linalg import splu
from .material_two_species import MaterialTwoSpecies


def local_geometry(coordinates, ids, xi, tau):
    shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
    derivative=np.stack((np.stack((-(1-tau),1-tau,tau,-tau),axis=-1),
                         np.stack((-(1-xi),-xi,xi,1-xi),axis=-1)),axis=-1)
    points=np.einsum('qn,cna->cqa',shape,coordinates[ids])
    jac=np.einsum('qnb,cna->cqab',derivative,coordinates[ids])
    gradient=np.einsum('qnb,cqba->cqna',derivative,np.linalg.inv(jac))
    return shape,points,gradient,np.linalg.det(jac),np.linalg.inv(jac)


class CompactMagneticField:
    def __init__(self,fluid,*,poloidal_tesla,toroidal_tesla):
        self.fluid=fluid;self.bp=float(poloidal_tesla);self.bt=float(toroidal_tesla)
        if not np.isfinite(self.bp+self.bt) or fluid.domain.radial_grading!=1:
            raise ValueError('Finite amplitudes and the documented uniform radial reference labels required')

    def __call__(self,xi,tau):
        f=self.fluid;d=f.domain
        shape,x,grad,det,inv=local_geometry(f.initial.coordinates,f.ids,xi,tau)
        s=(np.repeat(np.arange(d.nr),d.nt)[:,None]+xi[None,:])/d.nr
        h=1-s*s;r=x[...,0]
        ds=inv[...,0,:]/d.nr
        derivative=self.bp*(r[...,None]*h[...,None]**4*np.array([1.,0.])
            -4*r[...,None]**2*s[...,None]*h[...,None]**3*ds)
        result=np.zeros(r.shape+(3,),dtype=np.result_type(r,xi,tau))
        result[...,0]=-derivative[...,1]/r
        result[...,2]=derivative[...,0]/r
        result[...,1]=self.bt*r/d.radius*h**3
        return result


class MaterialMagnetic:
    def __init__(self,fluid,reference_field):
        self.fluid=fluid;self.reference_field=reference_field
        q,_=np.polynomial.legendre.leggauss(fluid.domain.quadrature)
        xi,tau=np.meshgrid(.5*(q+1),.5*(q+1),indexing='ij')
        self.xi=xi.ravel();self.tau=tau.ravel()
        self.initial_field=np.asarray(reference_field(self.xi,self.tau),float)
        if self.initial_field.shape!=fluid.particles.shape+(3,) or not np.isfinite(self.initial_field).all():
            raise ValueError('Finite reference magnetic field at material quadrature required')
        self.reference_volume=fluid.reference_geometry['volume']

    def response(self,state,fields=None):
        f=self.fluid;g=f.fields(state) if fields is None else fields
        deform=f.deformation(state);jac=np.linalg.det(deform)
        if np.any(jac<=0) or not np.isfinite(jac).all():raise ValueError('Positive full material deformation required')
        y=np.einsum('cqab,cqb->cqa',deform,self.initial_field);b=y/jac[...,None]
        scale=self.reference_volume/(MU0*jac);square=np.sum(y*y,axis=-1)
        reference_gradient=f.reference_geometry['gradient']
        directional=np.einsum('cqa,cqna->cqn',self.initial_field[...,[0,2]],reference_gradient)
        angle_gradient=np.einsum('cqna,cn->cqa',reference_gradient,state.angle[f.ids])
        angular_factor=np.sum(angle_gradient*self.initial_field[...,[0,2]],axis=-1)
        angular_factor+=self.initial_field[...,1]/f.reference_geometry['points'][...,0]
        local=np.einsum('cq,cqa,cqn->cna',scale,y[...,[0,2]],directional)
        local[...,0]+=np.einsum('cq,qn->cn',scale*y[...,1]*angular_factor,f.shape)
        local-=.5*np.einsum('cq,cqna->cna',scale*square,g['divergence'])
        force=np.zeros_like(state.coordinates);np.add.at(force,f.ids.ravel(),-local.reshape(-1,2));force[f.domain.axis,0]=0
        local_torque=-np.einsum('cq,cqn->cn',scale*y[...,1]*g['points'][...,0],directional)
        torque=np.zeros(f.nodes);np.add.at(torque,f.ids.ravel(),local_torque.ravel())
        return dict(magnetic_field=b,magnetic_energy=float(.5*np.sum(scale*square)),
            magnetic_force=force,magnetic_torque=torque)

    def _field_at(self,state,xi,tau):
        f=self.fluid
        shape,ref,grad0,det0,inv0=local_geometry(f.initial.coordinates,f.ids,xi,tau)
        _,x,grad,det,inv=local_geometry(state.coordinates,f.ids,xi,tau)
        drdz=np.einsum('cqnb,cna->cqab',grad0,state.coordinates[f.ids])
        twist=np.einsum('cqna,cn->cqa',grad0,state.angle[f.ids])
        b0=self.reference_field(xi,tau);r=x[...,0];radius=ref[...,0]
        jac=det*r/(det0*radius)
        result=np.zeros(b0.shape,dtype=np.result_type(b0,x))
        result[...,[0,2]]=np.einsum('cqab,cqb->cqa',drdz,b0[...,[0,2]])
        result[...,1]=r*(np.sum(twist*b0[...,[0,2]],axis=-1)+b0[...,1]/radius)
        return result/jac[...,None],x,inv

    def samples(self,state):
        b,x,inv=self._field_at(state,self.xi,self.tau)
        epsilon=1e-25
        dx=self._field_at(state,self.xi+1j*epsilon,self.tau)[0].imag/epsilon
        dy=self._field_at(state,self.xi,self.tau+1j*epsilon)[0].imag/epsilon
        gradient=np.einsum('cqab,cqbd->cqad',np.stack((dx,dy),axis=-1),inv)
        current=np.stack((-gradient[...,1,1],gradient[...,0,1]-gradient[...,2,0],
            gradient[...,1,0]+b[...,1]/x[...,0]),axis=-1)/MU0
        divergence=gradient[...,0,0]+b[...,0]/x[...,0]+gradient[...,2,1]
        # Weak Ampere reconstruction includes all internal face distributions.
        # Regular test vectors are (r*N,0,0), (0,r*N,0), (0,0,N).
        f=self.fluid;g=f.geometry(state.coordinates);r=x[...,0];gradient_n=g['gradient'];shape=f.shape
        loads=np.stack((r[...,None]*gradient_n[...,1]*b[...,1,None],
            -r[...,None]*gradient_n[...,1]*b[...,0,None]+(2*shape[None,:,:]+r[...,None]*gradient_n[...,0])*b[...,2,None],
            -gradient_n[...,0]*b[...,1,None]),axis=-1)
        local=np.einsum('cq,cqna->cna',g['volume']/MU0,loads)
        load=np.zeros((f.nodes,3));np.add.at(load,f.ids.ravel(),local.reshape(-1,3))
        radial_mass=splu(f.gram(g['volume']*r*r));axial_mass=splu(f.gram(g['volume']))
        coefficients=np.column_stack((radial_mass.solve(load[:,0]),radial_mass.solve(load[:,1]),axial_mass.solve(load[:,2])))
        projected=np.einsum('qn,cna->cqa',shape,coefficients[f.ids]);projected[...,:2]*=r[...,None]
        return dict(magnetic_field=b,current_density=projected,cell_current_density=current,
            current_coefficients=coefficients,magnetic_divergence=divergence)


class MaterialIdealPlasma:
    def __init__(self,fluid,magnetic):
        if magnetic.fluid is not fluid:raise ValueError('Common moving material fluid required')
        self.fluid=fluid;self.magnetic=magnetic;self.initial=fluid.initial

    def fields(self,state):
        fields=self.fluid.fields(state);magnetic=self.magnetic.response(state,fields)
        return dict(fields,**magnetic,fluid_energy=fields['energy'],total_energy=fields['energy']+magnetic['magnetic_energy'])

    def derivative(self,state):
        result=self.fluid.derivative(state);magnetic=self.magnetic.response(state)
        result.momentum+=magnetic['magnetic_force'];result.angular_momentum+=magnetic['magnetic_torque']
        return result

    def step(self,state,dt):
        return MaterialTwoSpecies.step(self,state,dt)
