"""Local finite-mass Lorentz/Ampere coupling of all reacting species moments."""
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from .species_field_exponential import advance as exponential_advance


class AtomicElectromagnetic:
    """Time-centred Maxwell source; spatial curls/Gauss transport are separate."""
    def __init__(self, source, reference_density, magnetic_field, *, method='midpoint'):
        self.source=source;self.density=float(reference_density)
        self.magnetic=np.asarray(magnetic_field,dtype=float).copy()
        if self.density<=0 or not np.isfinite(self.density) or self.magnetic.shape!=(3,) or not np.isfinite(self.magnetic).all():
            raise ValueError('Positive reference density and finite laboratory magnetic field required')
        self.magnetic.setflags(write=False)
        self.charge=np.r_[source.network.charge,-1.]*EV
        self.frame=np.asarray(getattr(source,'reference_velocity',np.zeros(3)))
        self.electric_scale=np.sqrt(EPS0/(self.density*EV))
        if method not in ('midpoint','exponential'):raise ValueError('Choose midpoint or exponential field source')
        self.method=method

    def electromagnetic_step(self,state,electric,dt):
        """Advance the fixed-population Lorentz/Ampere equations and their work."""
        electric=np.asarray(electric,dtype=float)
        if not np.isfinite(dt) or dt<=0 or electric.shape!=(3,) or not np.isfinite(electric).all():
            raise ValueError('Positive finite step and laboratory electric field required')
        s=self.source;f=s.fields(state);number=f['numbers']
        root=np.sqrt(s.mass/s.neutral_mass)*np.sqrt(number)
        h0=root[:,None]*(f['velocity']+self.frame)/s.speed_scale
        e0=self.electric_scale*electric
        omega=self.charge*np.sqrt(self.density/(EPS0*s.mass))*np.sqrt(number)
        if self.method=='exponential':
            active=(number>0)&(self.charge!=0);h1=h0.copy()
            h1[active],e1=exponential_advance(h0[active],e0,omega[active],
                (self.charge[active]/s.mass[active])[:,None]*self.magnetic,dt)
            change=h1-h0;average_h=.5*(h1+h0)
        else:
            change,e1,average_h=self.midpoint(h0,e0,omega,dt)
        result=state.copy();delta_p=root[:,None]*change
        result[s.momentum]+=delta_p
        # The integrated species work is its exact kinetic-energy difference,
        # transformed to the stored fixed frame; no global budget is imposed.
        result[s.energy]+=np.sum(average_h*change,axis=1)-delta_p@(self.frame/s.speed_scale)
        new_electric=e1/self.electric_scale
        s.fields(result)
        if not np.isfinite(new_electric).all():raise ValueError('Nonfinite electromagnetic state')
        return result,new_electric

    def midpoint(self,h0,e0,omega,dt):
        s=self.source
        half=dt/2;rotation=half*(self.charge/s.mass)[:,None]*self.magnetic
        # R=(I-C)^-1, C h=h cross rotation. This is a matrix inverse,
        # not velocity damping or a projection of an evolved state.
        cross=np.zeros((len(omega),3,3));x,y,z=rotation.T
        cross[:,0,1]=z;cross[:,0,2]=-y;cross[:,1,0]=-z
        cross[:,1,2]=x;cross[:,2,0]=y;cross[:,2,1]=-x
        inverse=(np.eye(3)+cross+rotation[:,:,None]*rotation[:,None,:])/(1+np.sum(rotation**2,axis=1))[:,None,None]
        rotated=np.einsum('sij,sj->si',inverse,h0)
        operator=np.eye(3)+half**2*np.einsum('s,sij->ij',omega**2,inverse)
        midpoint_e=np.linalg.solve(operator,e0-half*np.einsum('s,si->i',omega,rotated))
        midpoint_h=rotated+half*omega[:,None]*np.einsum('sij,j->si',inverse,midpoint_e)
        change=2*(midpoint_h-h0);e1=2*midpoint_e-e0
        return change,e1,midpoint_h

    def step(self,state,electric,time,dt,*,sources=True,heat=None):
        if not sources:return self.electromagnetic_step(state,electric,dt)
        def material(values,start,duration):
            a=0. if heat is None else heat(start)
            b=0. if heat is None else heat(start+duration/2)
            return self.source.step(values,duration,self.density,initial_heat=a,midpoint_heat=b)
        first=material(state,time,dt/2)
        coupled,field=self.electromagnetic_step(first,electric,dt)
        return material(coupled,time+dt/2,dt/2),field

    def observe(self,state,electric):
        s=self.source;f=s.fields(state);row=s.observe(state,self.density)
        electric=np.asarray(electric);field_energy=.5*EPS0*(electric@electric)/(self.density*EV)
        field_momentum=EPS0*np.cross(electric,self.magnetic)/(self.density*s.momentum_scale)
        row.update(electric_v_m=electric,electric_energy_ev=field_energy,
            total_energy_budget_ev=row['energy_budget_ev']+field_energy,
            total_momentum_scaled=row['momentum_scaled']+field_momentum,
            plasma_frequency_rad_s=float(np.sqrt(self.density/EPS0*np.sum(self.charge**2/s.mass*f['numbers']))),
            maximum_cyclotron_frequency_rad_s=float(np.max(np.where(f['numbers']>0,abs(self.charge/s.mass)*np.linalg.norm(self.magnetic),0.))))
        return row
