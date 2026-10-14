"""Simultaneous moving fluid and finite-mass toroidal magnetic/current dynamics."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,mu_0 as MU0
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from scipy.linalg import cho_factor,cho_solve


@dataclass
class MaterialCurrentState:
    fluid: object
    magnetic: np.ndarray

    def added(self,other,factor):return MaterialCurrentState(self.fluid.added(other.fluid,factor),self.magnetic+factor*other.magnetic)


def positive_solve(matrix,value):
    diagonal=np.sqrt(np.diag(matrix))
    if np.any(diagonal<=0) or not np.isfinite(diagonal).all():raise ValueError('Positive finite physical metric required')
    scaled=matrix/diagonal[:,None]/diagonal[None,:]
    return cho_solve(cho_factor(scaled,lower=True),value/diagonal)/diagonal


class MaterialCarrierCurrent:
    def __init__(self,fluid,space):
        self.fluid=fluid;self.space=space;self.mu=fluid.ion_mass*ME/fluid.mass;self.delta=(fluid.ion_mass-ME)/fluid.mass
        self.volume0=fluid.reference_geometry['volume'].ravel();self.particles=fluid.particles.ravel()
        self.shape=fluid.particles.shape
        ref=fluid.reference_geometry['points'].reshape(-1,2)
        if ref.shape[0]!=len(space.reference_r) or not np.allclose(ref,np.column_stack((space.reference_r,space.reference_z)),rtol=1e-13,atol=1e-16):
            raise ValueError('Identical fluid/field reference quadrature required')
        if not np.array_equal(fluid.gammas,np.array([5/3,5/3])):raise ValueError('This isentropic sector uses gamma=5/3 for both species')

    def geometry(self,state,fields=None):
        f=self.fluid;g=f.fields(state) if fields is None else fields
        if np.any(state.angle!=0) or np.any(state.angular_momentum!=0):raise ValueError('Toroidal-field invariant sector requires zero bulk rotation')
        ent=state.entropy.reshape(-1,2)
        tolerance=128*np.finfo(float).eps*max(1.,float(abs(ent).max()))
        if np.max(np.ptp(ent,axis=0))>tolerance or np.max(abs(ent[:,0]-ent[:,1]))>tolerance:
            raise ValueError('Constant equal species entropies required; nonzero species transport/exchange is not omitted')
        deformation=np.einsum('cqnb,cna->cqab',f.reference_geometry['gradient'],state.coordinates[f.ids]).reshape(-1,2,2)
        inverse=np.linalg.inv(deformation);volume=g['volume'].ravel();r=g['points'][...,0].ravel();ratio=volume/self.volume0
        er=diags(inverse[:,0,0])@self.space.er+diags(inverse[:,1,0])@self.space.ez
        ez=diags(inverse[:,0,1])@self.space.er+diags(inverse[:,1,1])@self.space.ez
        factor=r/self.space.reference_r/ratio
        bmap=diags(factor)@self.space.magnetic;curl=diags(factor)@self.space.curl
        mass=(er.T@diags(self.particles)@er+ez.T@diags(self.particles)@ez).tocsc();solve=splu(mass).solve
        h=(bmap.T@diags(volume/MU0)@bmap).toarray();d=-(bmap.T@diags(volume/MU0)@curl).tocsr()
        gradient=np.einsum('cqnb,cna->cqab',g['gradient'],g['velocity'][f.ids]).reshape(-1,2,2)
        radial_rate=g['sampled_velocity'][...,0].ravel()/r;divergence=np.trace(gradient,axis1=1,axis2=2)+radial_rate
        metric_rate=2*radial_rate-divergence
        hdot=(bmap.T@diags(volume*metric_rate/MU0)@bmap).toarray()
        ddot=-(bmap.T@diags(volume*metric_rate/MU0)@curl).tocsr()
        grdot=-diags(gradient[:,0,0])@er-diags(gradient[:,1,0])@ez
        gzdot=-diags(gradient[:,0,1])@er-diags(gradient[:,1,1])@ez
        return dict(fields=g,er=er,ez=ez,bmap=bmap,curl=curl,mass=mass,solve=solve,h=h,d=d,hdot=hdot,ddot=ddot,
            grdot=grdot,gzdot=gzdot,gradient=gradient,volume=volume,r=r)

    def project_initial(self,fluid_state,magnetic_field,current_density):
        """Initial physical B/drift least squares, never an evolved-state repair."""
        g=self.geometry(fluid_state);b=np.asarray(magnetic_field,float).reshape(-1)
        j=np.asarray(current_density,float).reshape(-1,2)
        if b.shape!=self.volume0.shape or j.shape!=(len(b),2) or not np.isfinite(b).all() or not np.isfinite(j).all():
            raise ValueError('Consistent finite analytic initial Btheta and poloidal current samples required')
        magnetic_load=g['bmap'].T@(g['volume']/MU0*b)
        # integral dN G.w_exact = integral dV G.J_exact/e, with no J/n divide.
        drift_load=g['er'].T@(g['volume']/E*j[:,0])+g['ez'].T@(g['volume']/E*j[:,1])
        operator=g['h']+self.mu/E**2*(g['d']@g['solve'](g['d'].T.toarray()))
        rhs=magnetic_load-self.mu/E*(g['d']@g['solve'](drift_load))
        return MaterialCurrentState(fluid_state,positive_solve(operator,np.asarray(rhs)))

    def response(self,state,*,rate=True):
        f=self.fluid;g=self.geometry(state.fluid);b=state.magnetic
        if b.shape!=(self.space.magnetic_count,) or not np.isfinite(b).all():raise ValueError('Finite complete magnetic coordinates required')
        h=g['h'];d=g['d'];w=g['solve'](-d.T@b/E)
        wr=g['er']@w;wz=g['ez']@w;physical_b=g['bmap']@b
        energy_b=float(.5*b@h@b);energy_j=float(.5*self.mu*w@(g['mass']@w))
        result=dict(g,drift_coefficients=w,drift=np.column_stack((wr,wz)),magnetic_field=physical_b,
            magnetic_energy=energy_b,current_energy=energy_j,total_energy=g['fields']['energy']+energy_b+energy_j)
        if not rate:return result
        curlw=g['curl']@w;omega=g['gradient'][:,0,1]-g['gradient'][:,1,0]
        coefficient=self.delta*(physical_b+self.mu/E*curlw)-self.mu/E*omega
        forcing=g['er'].T@(-self.particles*wz*coefficient)+g['ez'].T@(self.particles*wr*coefficient)
        mwdot=g['grdot'].T@(self.particles*wr)+g['gzdot'].T@(self.particles*wz)
        mwdot+=g['er'].T@(self.particles*(g['grdot']@w))+g['ez'].T@(self.particles*(g['gzdot']@w))
        rhs=forcing-self.mu/E**2*(g['ddot'].T@b)-self.mu/E*mwdot
        inverse_d=g['solve'](d.T.toarray());operator=h+self.mu/E**2*(d@inverse_d)
        bdot=positive_solve(operator,np.asarray(d@g['solve'](rhs)))
        electric=g['solve'](rhs-self.mu/E**2*(d.T@bdot))
        # Refine the original mixed equations, not a state or energy budget.
        # Recovering E by subtracting two large Schur terms can lose Faraday
        # accuracy at dilute carrier edges even when the Schur solve is sound.
        for refinement in range(4):
            r1=d@electric-h@bdot
            r2=rhs-self.mu/E**2*(d.T@bdot)-g['mass']@electric
            scale1=max(np.linalg.norm(d@electric)+np.linalg.norm(h@bdot),np.finfo(float).tiny)
            scale2=max(np.linalg.norm(rhs)+np.linalg.norm(self.mu/E**2*(d.T@bdot))+np.linalg.norm(g['mass']@electric),np.finfo(float).tiny)
            if np.linalg.norm(r1)<2e-13*scale1 and np.linalg.norm(r2)<2e-13*scale2:break
            correction=positive_solve(operator,r1+d@g['solve'](r2))
            electric+=g['solve'](r2-self.mu/E**2*(d.T@correction));bdot+=correction
        wdot=g['solve'](-(d.T@bdot+g['ddot'].T@b)/E-mwdot)
        faraday=h@bdot-d@electric
        ohm=g['mass']@electric-forcing-self.mu/E*(g['mass']@wdot)
        norm=lambda residual,*terms:float(np.linalg.norm(residual)/max(sum(np.linalg.norm(v) for v in terms),np.finfo(float).tiny))
        errors=dict(faraday=norm(faraday,h@bdot,d@electric),ohm=norm(ohm,g['mass']@electric,forcing,self.mu/E*(g['mass']@wdot)))
        if max(errors.values())>1e-10 or not np.isfinite(bdot).all():raise RuntimeError(f'Physical moving induction residual failed: {errors}')
        fields=g['fields'];vol=fields['volume'];bv=physical_b.reshape(self.shape);wv=np.stack((wr,wz),axis=-1).reshape(self.shape+(2,))
        local=np.einsum('cq,cqna->cna',vol*bv*bv/(2*MU0),fields['divergence'])
        local[...,0]-=np.einsum('cq,qn->cn',vol*bv*bv/(MU0*fields['points'][...,0]),f.shape)
        directional=np.einsum('cqa,cqna->cqn',wv,fields['gradient'])
        local+=self.mu*np.einsum('cq,cqa,cqn->cna',f.particles,wv,directional)
        force=np.zeros_like(state.fluid.coordinates);np.add.at(force,f.ids.ravel(),local.reshape(-1,2));force[f.domain.axis,0]=0
        current_power=self.mu*w@(g['mass']@wdot)+.5*self.mu*w@mwdot
        magnetic_power=b@h@bdot+.5*b@g['hdot']@b
        mechanical_power=float(np.sum(force*fields['velocity']))
        charge=self.space.gradient.T@(g['mass']@w)
        result.update(magnetic_rate=bdot,electric_reduced=electric,drift_rate=wdot,force=force,residuals=errors,
            magnetic_power=float(magnetic_power),current_power=float(current_power),mechanical_power=mechanical_power,
            pair_power_residual=float(current_power+magnetic_power+mechanical_power),
            weak_charge_residual=float(np.linalg.norm(charge)/max(np.linalg.norm(self.space.gradient)*np.linalg.norm(g['mass']@w),np.finfo(float).tiny)))
        return result

    def derivative(self,state):
        response=self.response(state);fluid=self.fluid.derivative(state.fluid);fluid.momentum+=response['force']
        return MaterialCurrentState(fluid,response['magnetic_rate'])

    def step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite step required')
        k1=self.derivative(state);k2=self.derivative(state.added(k1,.5*dt));k3=self.derivative(state.added(k2,.5*dt));k4=self.derivative(state.added(k3,dt))
        result=state.added(k1,dt/6).added(k2,dt/3).added(k3,dt/3).added(k4,dt/6)
        self.response(result,rate=False);return result
