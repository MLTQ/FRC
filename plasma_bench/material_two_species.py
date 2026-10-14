"""Two-species material fluid with vanishing density at a free vacuum interface."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


@dataclass
class MaterialSpeciesState:
    coordinates: np.ndarray
    momentum: np.ndarray
    angular_momentum: np.ndarray
    angle: np.ndarray
    entropy: np.ndarray

    def added(self,other,factor):
        return MaterialSpeciesState(*(getattr(self,k)+factor*getattr(other,k) for k in self.__dataclass_fields__))


class MaterialTwoSpecies:
    """General meridional Q1 motion, consistent inertia and two physical entropies."""
    def __init__(self,domain,coordinates,number_density,temperature_ev,ion_mass,*,gammas=(5/3,5/3),velocity=None,angular_velocity=None):
        self.domain=domain;self.ids=domain.cells[:domain.fluid_cell_count];self.shape=domain.shape;self.reference_gradient=domain.derivative
        self.ion_mass=float(ion_mass);self.mass=self.ion_mass+ME;self.gammas=np.asarray(gammas,float)
        if self.ion_mass<=ME or self.gammas.shape!=(2,) or np.any(self.gammas<=1):raise ValueError('Physical ion mass and two gamma>1 values required')
        self.active_r=np.flatnonzero(~domain.axis);self.nodes=domain.fluid_node_count
        coordinates=np.asarray(coordinates,float);geo=self.geometry(coordinates);self.reference_geometry=geo;n=np.broadcast_to(number_density,geo['volume'].shape).copy()
        temperature=np.broadcast_to(temperature_ev,n.shape+(2,)).copy()
        if not np.isfinite(n).all() or not np.isfinite(temperature).all() or np.any(n<=0) or np.any(temperature<=0):
            raise ValueError('Positive physical density and temperature at interior material quadrature required; no floors')
        self.particles=n*geo['volume'];self.particles.setflags(write=False)
        self.mass_matrix=self.gram(self.mass*self.particles)
        self.mass_factors=(splu(self.mass_matrix[self.active_r][:,self.active_r].tocsc()),splu(self.mass_matrix))
        v=np.zeros_like(coordinates) if velocity is None else np.asarray(velocity,float)
        omega=np.zeros(self.nodes) if angular_velocity is None else np.asarray(angular_velocity,float)
        if v.shape!=coordinates.shape or omega.shape!=(self.nodes,) or np.any(v[domain.axis,0]!=0) or not np.isfinite(v).all() or not np.isfinite(omega).all():
            raise ValueError('Finite physical meridional and angular velocities required')
        momentum=np.column_stack((self.mass_matrix@v[:,0],self.mass_matrix@v[:,1]));momentum[domain.axis,0]=0
        angular=self.gram(self.mass*self.particles*geo['points'][...,0]**2)@omega
        # Entropy convention T[eV]=exp(s)*(n/nref)^(gamma-1), nref=1/m3.
        entropy=np.log(temperature)-(self.gammas-1)*np.log(n)[...,None]
        self.initial=MaterialSpeciesState(coordinates.copy(),momentum,angular,np.zeros(self.nodes),entropy)

    def gram(self,weight):
        local=np.einsum('cq,qi,qj->cij',weight,self.shape,self.shape)
        return coo_matrix((local.ravel(),(np.repeat(self.ids,4,axis=1).ravel(),np.tile(self.ids,(1,4)).ravel())),shape=(self.nodes,self.nodes)).tocsc()

    def geometry(self,coordinates):
        d=self.domain;x=np.asarray(coordinates,float)
        if x.shape!=d.fluid_reference.shape or not np.isfinite(x).all() or np.any(x[d.axis,0]!=0) or np.any(x[~d.axis,0]<=0):
            raise ValueError('Finite conforming material geometry with exact regular axis required')
        corners=np.array([[[-1.,-1.],[1.,0.],[0.,0.],[0.,1.]],
            [[-1.,0.],[1.,-1.],[0.,1.],[0.,0.]],[[0.,0.],[0.,-1.],[1.,1.],[-1.,0.]],
            [[0.,-1.],[0.,0.],[1.,0.],[-1.,1.]]])
        det=np.linalg.det(np.einsum('qnb,cna->cqab',corners,x[self.ids]))
        allowed=np.zeros_like(det,dtype=bool);allowed[:d.nt,[0,3]]=True
        if np.any(det<0) or np.any((det==0)&~allowed):raise ValueError('Invalid material cell corner; no mesh repair')
        points=np.einsum('qn,cna->cqa',self.shape,x[self.ids])
        jac=np.einsum('qnb,cna->cqab',self.reference_gradient,x[self.ids]);det=np.linalg.det(jac)
        if np.any(det<=0) or np.any(points[...,0]<=0) or not np.isfinite(det).all():raise ValueError('Invalid interior material Jacobian/radius')
        gradient=np.einsum('qnb,cqba->cqna',self.reference_gradient,np.linalg.inv(jac))
        weight=d.qweight[None,:]*det;volume=2*np.pi*points[...,0]*weight
        divergence=gradient.copy();divergence[...,0]+=self.shape[None,:,:]/points[...,0,None]
        return dict(points=points,volume=volume,divergence=divergence,gradient=gradient,jacobian=det)

    def fields(self,state):
        if state.momentum.shape!=(self.nodes,2) or state.angular_momentum.shape!=(self.nodes,) or state.angle.shape!=(self.nodes,) or state.entropy.shape!=self.particles.shape+(2,):
            raise ValueError('Complete material canonical state required')
        if not all(np.isfinite(getattr(state,k)).all() for k in state.__dataclass_fields__) or np.any(state.momentum[self.domain.axis,0]!=0):
            raise ValueError('Finite canonical state and inactive radial-axis momentum required')
        g=self.geometry(state.coordinates);velocity=np.zeros_like(state.momentum)
        velocity[self.active_r,0]=self.mass_factors[0].solve(state.momentum[self.active_r,0]);velocity[:,1]=self.mass_factors[1].solve(state.momentum[:,1])
        radial=g['points'][...,0]
        if np.any(state.angular_momentum):
            angular_mass=self.gram(self.mass*self.particles*radial**2);omega=splu(angular_mass).solve(state.angular_momentum)
        else:omega=np.zeros(self.nodes)
        sampled_velocity=np.einsum('qn,cna->cqa',self.shape,velocity[self.ids]);sampled_omega=np.einsum('qn,cn->cq',self.shape,omega[self.ids])
        n=self.particles/g['volume'];temperature=np.exp(state.entropy+(self.gammas-1)*np.log(n)[...,None])
        if np.any(n<=0) or np.any(temperature<=0) or not np.isfinite(temperature).all():raise ValueError('Nonpositive/nonfinite material species, no clipping')
        specific=E*temperature/(self.gammas-1);pressure=n[...,None]*E*temperature
        internal=np.sum(self.particles[...,None]*specific,axis=(0,1))
        kinetic=.5*np.sum(state.momentum*velocity)+.5*state.angular_momentum@omega
        return dict(g,n=n,temperature=temperature,specific_internal=specific,pressure=pressure,
            velocity=velocity,sampled_velocity=sampled_velocity,omega=omega,sampled_omega=sampled_omega,
            internal_energy=internal,kinetic_energy=float(kinetic),energy=float(kinetic+internal.sum()))

    def deformation(self,state):
        """Full axisymmetric 3-D material gradient in local orthonormal bases."""
        g0=self.reference_geometry;g=self.geometry(state.coordinates)
        drdz=np.einsum('cqnb,cna->cqab',g0['gradient'],state.coordinates[self.ids])
        twist=np.einsum('cqna,cn->cqa',g0['gradient'],state.angle[self.ids])
        result=np.zeros(self.particles.shape+(3,3));r=g['points'][...,0]
        result[...,0,0]=drdz[...,0,0];result[...,0,2]=drdz[...,0,1]
        result[...,2,0]=drdz[...,1,0];result[...,2,2]=drdz[...,1,1]
        result[...,1,0]=r*twist[...,0];result[...,1,2]=r*twist[...,1]
        result[...,1,1]=r/g0['points'][...,0]
        return result

    def derivative(self,state):
        f=self.fields(state)
        # -dU/dx = integral p div(delta x) dV; pressure vanishes at the
        # material vacuum edge through fixed compact carrier/entropy profiles.
        local=np.einsum('cq,cqna->cna',f['volume']*np.sum(f['pressure'],axis=-1),f['divergence'])
        # Canonical angular momentum is fixed: -dKtheta/dx gives centrifugal force.
        local[...,0]+=np.einsum('cq,qn->cn',self.mass*self.particles*f['points'][...,0]*f['sampled_omega']**2,self.shape)
        force=np.zeros_like(state.coordinates);np.add.at(force,self.ids.ravel(),local.reshape(-1,2));force[self.domain.axis,0]=0
        return MaterialSpeciesState(f['velocity'],force,np.zeros_like(state.angular_momentum),f['omega'],np.zeros_like(state.entropy))

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite material step required')
        k1=self.derivative(state);k2=self.derivative(state.added(k1,.5*dt));k3=self.derivative(state.added(k2,.5*dt));k4=self.derivative(state.added(k3,dt))
        result=MaterialSpeciesState(*(getattr(state,k)+dt/6*(getattr(k1,k)+2*getattr(k2,k)+2*getattr(k3,k)+getattr(k4,k)) for k in state.__dataclass_fields__))
        self.fields(result);return result
