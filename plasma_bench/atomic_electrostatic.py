"""Charge-continuous longitudinal fields and work for spatial species moments."""
import numpy as np
from scipy.constants import elementary_charge as EV, epsilon_0 as EPS0
from .atomic_transport import AtomicTransport
from .species_euler_flux import hllc


class AtomicElectrostatic(AtomicTransport):
    """Reflecting interval, Gauss field, and matched numerical charge-flux work."""
    def __init__(self,*args,left_electric=0.,maximum_plasma_phase=.02,**kwargs):
        super().__init__(*args,**kwargs)
        if not np.isfinite([left_electric,maximum_plasma_phase]).all() or maximum_plasma_phase<=0:
            raise ValueError('Finite boundary field and positive phase step required')
        self.left_electric=float(left_electric);self.maximum_plasma_phase=float(maximum_plasma_phase)
        self.charge=np.r_[self.source.network.charge,-1.]

    def electric(self,state):
        charge=EV*self.density*(state[:,self.number]@self.charge)
        return np.r_[self.left_electric,self.left_electric+self.dx/EPS0*np.cumsum(charge)]

    def timestep(self,state):
        acoustic=super().timestep(state)
        omega=np.sqrt(EV**2*self.density/EPS0*np.max(state[:,self.number]@(self.charge**2/self.source.mass)))
        return min(acoustic,self.maximum_plasma_phase/omega) if omega else acoustic

    def face_flux(self,state):
        p=self.primitive(state)
        def reflect(values):
            ghost=values.copy();ghost[...,4]=-values[...,4]-2*self.frame[2];return ghost
        extended=np.concatenate((reflect(p[:1]),p,reflect(p[-1:])))
        backward=extended[1:-1]-extended[:-2];forward=extended[2:]-extended[1:-1]
        slope=np.where(backward*forward>0,np.sign(backward)*np.minimum(abs(backward),abs(forward)),0.)
        lo,hi=p-.5*slope,p+.5*slope
        left=np.concatenate((reflect(lo[:1]),hi));right=np.concatenate((lo,reflect(hi[-1:])))
        flux=hllc(left,right,self.mass_ratio,self.frame[2],gamma=self.gamma)*self.source.speed_scale
        # Exact impermeable stationary-wall condition. Only normal pressure
        # momentum flux remains; transform its zero LAB energy flux to the frame.
        flux[[0,-1],:,0:3]=0.
        flux[[0,-1],:,4]=-self.frame[2]*flux[[0,-1],:,3]
        return flux

    def derivative(self,state):
        s=self.source;flux=self.face_flux(state);field=self.electric(state)
        rate=-(flux[1:]-flux[:-1])/self.dx;result=np.zeros_like(state)
        result[:,self.number]=rate[...,0];result[:,s.momentum]=rate[...,1:4];result[:,s.energy]=rate[...,4]
        centre=.5*(field[:-1]+field[1:]);numbers=state[:,self.number]
        force=(EV/s.momentum_scale)*self.charge*numbers*centre[:,None]
        result[:,s.momentum[:,2]]+=force
        work=.5*self.charge*(field[:-1,None]*flux[:-1,:,0]+field[1:,None]*flux[1:,:,0])
        result[:,s.energy]+=work-(s.momentum_scale/EV)*self.frame[2]*s.speed_scale*force
        current=EV*self.density*(flux[:,:,0]@self.charge)
        field_rate=np.r_[0.,self.dx*EV*self.density/EPS0*np.cumsum(result[:,self.number]@self.charge)]
        field_power=EPS0/(EV*self.density)*np.trapezoid(field*field_rate,dx=self.dx)
        self.last_diagnostic=dict(
            ampere_residual_v_m_s=float(abs(field_rate+current/EPS0).max()),
            electric_work_residual_ev_m_s=float(self.dx*work.sum()+field_power))
        return result

    def totals(self,state):
        row=super().totals(state);s=self.source;field=self.electric(state)
        lab=row['energy_ev']+(s.momentum_scale/EV)*s.speed_scale*(self.frame@row['momentum_scaled'])
        lab+=.5*(s.neutral_mass/EV)*s.speed_scale**2*(self.frame@self.frame)*row['mass_over_neutral']
        field_energy=EPS0/(2*EV*self.density)*np.trapezoid(field**2,dx=self.dx)
        exact_linear_energy=EPS0*self.dx/(6*EV*self.density)*np.sum(field[:-1]**2+field[:-1]*field[1:]+field[1:]**2)
        gauss=EPS0*np.diff(field)/self.dx-EV*self.density*(state[:,self.number]@self.charge)
        row.update(particle_lab_energy_ev=lab,electric_energy_ev=field_energy,
            total_energy_ev=lab+field_energy,linear_field_total_energy_ev=lab+exact_linear_energy,
            maximum_gauss_residual_over_reference_charge=float(abs(gauss).max()/(EV*self.density)))
        return row
