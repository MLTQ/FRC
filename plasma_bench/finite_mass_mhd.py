"""Explicit unsplit bounded finite-mass reference evolution with atomic steps."""
import numpy as np
from .finite_mass_state import FiniteMassState,evaluate
from .finite_mass_rhs import rhs


class FiniteMassMHD:
    """Research reference: resolve stiff current modes with an explicit time cap."""
    def __init__(self,mesh,initial,parameters,*,max_timestep,source=None,initial_time=0.):
        if not np.isfinite(max_timestep) or max_timestep<=0:raise ValueError('Explicit positive maximum time step required')
        if not np.isfinite(initial_time):raise ValueError('Finite initial time required')
        evaluate(mesh,initial,parameters)
        self.mesh,self.parameters=mesh,parameters;self.state=initial.copy()
        self.max_timestep,self.source=max_timestep,source
        self.time=float(initial_time);self.steps=0;self.energy_input=0.;self.max_linear_residual=0.;self.linear_iterations=0
        self.initial_integrals=self.integrals()

    def derivative(self,state,time):return rhs(self.mesh,state,self.parameters,time=time,source=self.source)

    def step(self,dt=None):
        dt=self.max_timestep if dt is None else dt
        if not np.isfinite(dt) or dt<=0 or dt>self.max_timestep*(1+1e-12):
            raise ValueError('Positive step within the explicit requested time cap required')
        old=self.state;t=self.time
        a=self.derivative(old,t);b=self.derivative(old.add(a.state,.5*dt),t+.5*dt)
        c=self.derivative(old.add(b.state,.5*dt),t+.5*dt);d=self.derivative(old.add(c.state,dt),t+dt)
        stages=(a,b,c,d);weights=(1/6,1/3,1/3,1/6)
        delta=FiniteMassState(*(sum(w*s.state.arrays()[k] for w,s in zip(weights,stages)) for k in range(4)))
        trial=old.add(delta,dt);evaluate(self.mesh,trial,self.parameters)
        def power(stage):
            boundary=np.sum(stage.radial_flux[-1,:,3]-stage.radial_flux[0,:,3])
            boundary+=np.sum(stage.axial_flux[:,-1,3]-stage.axial_flux[:,0,3])
            return float(np.sum(self.mesh.volumes*stage.fluid_source[...,3])-boundary)
        work=dt*sum(w*power(s) for w,s in zip(weights,stages))
        # Physical state, time and every accepted budget commit together.
        self.state=trial;self.time+=dt;self.steps+=1;self.energy_input+=work
        self.max_linear_residual=max(self.max_linear_residual,*(s.linear_residual for s in stages))
        self.linear_iterations+=sum(s.linear_iterations for s in stages)
        return dt

    def advance(self,end_time,*,max_steps=100000):
        if not np.isfinite(end_time) or end_time<self.time:raise ValueError('Finite nondecreasing target time required')
        count=0
        while self.time<end_time:
            remaining=end_time-self.time
            if remaining<=4*np.finfo(float).eps*max(abs(end_time),1e-300):self.time=end_time;break
            if count>=max_steps:raise RuntimeError('Finite-mass step budget exceeded; accepted partial state retained')
            self.step(min(remaining,self.max_timestep));count+=1
        return self.state

    def integrals(self):
        return {name:float(np.sum(self.mesh.volumes*self.state.conserved[...,column]))
                for name,column in (('mass_kg',0),('axial_momentum_kg_m_s',2),('energy_j',3),('angular_momentum_kg_m2_s',4))}

    def diagnostics(self):
        result=self.integrals();ev=evaluate(self.mesh,self.state,self.parameters)
        from scipy.constants import elementary_charge as E
        n=ev.model.electron.density;u=self.state.electron_entropy
        ion_n=self.state.conserved[...,0]/self.parameters.masses.total
        pi=ev.components['ion']*(self.parameters.gamma_i-1)
        ion_w=pi/(ion_n*E)*(ion_n/self.parameters.density_scale)**(1-self.parameters.gamma_i)
        result.update(time_s=self.time,steps=self.steps,external_energy_input_j=self.energy_input,
            energy_budget_residual_j=result['energy_j']-self.initial_integrals['energy_j']-self.energy_input,
            min_electron_pressure_pa=float(ev.model.electron.pressure(self.state.electron_entropy).min()),
            min_ion_pressure_pa=float((ev.components['ion']*(self.parameters.gamma_i-1)).min()),
            min_number_density_m3=float(ev.model.electron.density.min()),
            max_linear_residual=self.max_linear_residual,linear_iterations=self.linear_iterations,
            electron_entropy_proxy_j=float(np.sum(ev.model.mapping.volumes*u)),
            electron_nodal_log_entropy_particles=float(np.sum(ev.model.mapping.volumes*n*np.log(u/(n*E)))),
            ion_cell_log_entropy_particles=float(np.sum(self.mesh.volumes*ion_n*np.log(ion_w))),
            scope='Explicit bounded positive-density finite-mass reference; no vacuum/transport/formation validation')
        return result
