"""Independent physical response and energy derivatives of saved full rates."""
import numpy as np
from scipy.constants import elementary_charge as E,mu_0 as MU0
from .trajectory_tangent_audit import expand,components
from .trajectory_density_prefix import fingerprint


def response(sample,states,rates,amplitude):
    minus=expand(sample,states['minus_state'],rates['minus']['fluid_rate'],rates['minus']['magnetic_rate'])
    plus=expand(sample,states['plus_state'],rates['plus']['fluid_rate'],rates['plus']['magnetic_rate'])
    fields=(plus[0]-minus[0])/(2*amplitude);pressure=(plus[1]-minus[1])/(2*amplitude)
    current=sample.electric((rates['plus']['current_rate']-rates['minus']['current_rate'])/(2*amplitude))
    electric=sample.electric((rates['plus']['electric']-rates['minus']['electric'])/(2*amplitude))
    return components(fields,pressure,current,electric),{k:fingerprint(v) for k,v in
        dict(rhs_fields=fields,rhs_pressure=pressure,current_rate=current,electric=electric).items()}


def power(sample,state,magnetic,current,rates,bias,ion_mass,electron_mass):
    n,v,_=sample.fluid(state)
    fields,pt,_=expand(sample,state,rates['fluid_rate'],rates['magnetic_rate'])
    nt=fields[...,0];vt=fields[...,1:4];bt=fields[...,6:9]
    j=sample.electric(current);jt=sample.electric(rates['current_rate'])
    b=sample.magnetic(magnetic);b[...,2]+=bias
    mass=ion_mass+electron_mass;reduced=ion_mass*electron_mass/mass
    integrate=lambda value:float(np.sum(sample.weight*value))
    powers=dict(fluid=integrate(mass*(.5*nt*np.sum(v*v,axis=-1)+n*np.sum(v*vt,axis=-1))+np.sum(pt/(sample.d.thermal.gammas-1),axis=-1)),
                current=integrate(reduced/E**2*(np.sum(j*jt,axis=-1)/n-.5*np.sum(j*j,axis=-1)*nt/n**2)),
                magnetic=integrate(np.sum(b*bt,axis=-1)/MU0))
    total=sum(powers.values());scale=sum(abs(x) for x in powers.values())
    return dict(powers_W=powers,total_power_W=total,absolute_exchange_scale_W=scale,
        total_relative_to_exchange=abs(total)/scale if scale else None,
        particle_rate=integrate(nt),axial_momentum_rate=integrate(mass*(nt*v[...,2]+n*vt[...,2])),
        angular_momentum_rate=integrate(mass*sample.R*(nt*v[...,1]+n*vt[...,1])),
        entropy_density_rates=[integrate(sample.scalar(rates['fluid_rate'][...,k])) for k in (4,5)],
        input_current_fingerprint=fingerprint(j))
