"""Experimental ideal material plasma evolving against a massless infinite vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .waveforms import DC
from .moving_field import MovingField
from .moving_advection import MaterialFlux


@dataclass
class MovingState:
    coordinates: np.ndarray
    velocity: np.ndarray

    def copy(self):
        return MovingState(self.coordinates.copy(),self.velocity.copy())


class MovingPlasma:
    """Adiabatic ideal reference; no remapping, shock capture or transport closure."""
    def __init__(self, domain, initial_a, *, density=1., pressure=1., velocity=None,
                 bias=0., gamma=5/3, cfl=.15, max_timestep=None):
        if not np.isfinite(gamma) or gamma <= 1 or not np.isfinite(cfl) or not 0 < cfl <= .5:
            raise ValueError("Require gamma>1 and 0<CFL<=0.5")
        if max_timestep is not None and (not np.isfinite(max_timestep) or max_timestep <= 0):
            raise ValueError("Source/accuracy time cap must be finite and positive")
        self.domain, self.gamma, self.cfl, self.max_timestep = domain,gamma,cfl,max_timestep
        self.bias = DC(float(bias)) if np.isscalar(bias) else bias
        if not isinstance(self.bias,DC) and max_timestep is None:
            raise ValueError("A time-dependent prescribed field requires an explicit source time cap")
        self.field = MovingField(domain)
        coordinates = domain.fluid_reference.copy()
        a = initial_a(coordinates[:,0],coordinates[:,1]) if callable(initial_a) else initial_a
        self.flux = MaterialFlux(domain,np.broadcast_to(a,(domain.fluid_node_count,)))
        self.state = MovingState(coordinates,np.zeros_like(coordinates) if velocity is None else np.asarray(velocity,float).copy())
        if self.state.velocity.shape != coordinates.shape or not np.all(np.isfinite(self.state.velocity)) or np.any(self.state.velocity[domain.axis,0] != 0):
            raise ValueError("Finite nodal velocity and zero radial axis velocity required")
        batch = domain.elements(coordinates)[0]; count = domain.fluid_cell_count
        r,z = batch.points[:count,:,0],batch.points[:count,:,1]
        rho = np.broadcast_to(density(r,z) if callable(density) else density,r.shape)
        p = np.broadcast_to(pressure(r,z) if callable(pressure) else pressure,r.shape)
        if np.any(rho <= 0) or np.any(p < 0) or not np.all(np.isfinite(rho+p)):
            raise ValueError("Require positive material density and nonnegative pressure; no floors")
        weight = 2*np.pi*batch.weights[:count]*r
        volume = weight.sum(axis=1)
        self.cell_mass = np.sum(weight*rho,axis=1)
        self.pressure_invariant = np.sum(weight*p,axis=1)*volume**(gamma-1)
        local_mass = np.einsum("cq,cqn->cn",weight*rho,batch.geometry_shape[:count])
        self.node_mass = np.zeros(domain.fluid_node_count)
        np.add.at(self.node_mass,domain.cells[:count].ravel(),local_mass.ravel())
        if np.any(self.node_mass <= 0):
            raise ValueError("Every material velocity node must carry positive physical mass")
        self.time, self.steps, self.source_work = 0.,0,0.
        self.initial = self.diagnostics()

    def evaluate(self, state=None, time=None):
        state = self.state if state is None else state
        time = self.time if time is None else time
        d = self.domain
        bias,rate = float(self.bias.value(time)),float(self.bias.derivative(time))
        if not np.isfinite(bias+rate):
            raise ValueError("Prescribed bias and derivative must be finite")
        a = self.flux.coefficients(state.coordinates)
        magnetic = self.field.solve(state.coordinates,a-.5*bias)
        volume = d.volumes(magnetic.elements)
        pressure = self.pressure_invariant/volume**self.gamma
        internal = float(np.sum(self.pressure_invariant*volume**(1-self.gamma)/(self.gamma-1)))
        force = -magnetic.geometry_gradient.copy()
        force[:,0] -= self.flux.pullback(state.coordinates,magnetic.load[:d.fluid_node_count])
        pressure_force = pressure[:,None,None]*d.volume_gradients(magnetic.elements)
        np.add.at(force,d.cells[:d.fluid_cell_count].ravel(),pressure_force.reshape(-1,2))
        force[d.axis,0] = 0.  # Constrained symmetry reaction performs no work.
        source_power = -.5*rate*float(magnetic.load[:d.fluid_node_count].sum())
        return {"magnetic":magnetic,"volume":volume,"pressure":pressure,"internal_energy":internal,
                "force":force,"source_power":source_power,"bias":bias}

    def timestep_limit(self, evaluation=None):
        e = self.evaluate() if evaluation is None else evaluation
        d = self.domain;count = d.fluid_cell_count;batch = e["magnetic"].elements[0]
        br,bz = self.field.fields(batch,e["magnetic"].induced,e["bias"])
        rho = self.cell_mass/e["volume"]
        speed = np.sqrt(self.gamma*e["pressure"]/rho+np.max(br[:count]**2+bz[:count]**2,axis=1)/(MU0*rho))
        nodes = self.state.coordinates[d.cells[:count]]
        perimeter = np.linalg.norm(np.roll(nodes,-1,axis=1)-nodes,axis=-1).sum(axis=1)
        width = 2*batch.weights[:count].sum(axis=1)/perimeter
        velocity = self.state.velocity[d.cells[:count]]
        relative = np.max(np.linalg.norm(velocity-velocity.mean(axis=1)[:,None,:],axis=-1),axis=1)
        inverse = np.max((speed+relative)/width)
        limit = np.inf if inverse == 0 else self.cfl/inverse
        return float(limit if self.max_timestep is None else min(limit,self.max_timestep))

    def step(self, dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Step must be finite and positive")
        old = self.evaluate()
        half_velocity = self.state.velocity+.5*dt*old["force"]/self.node_mass[:,None]
        proposal = MovingState(self.state.coordinates+dt*half_velocity,half_velocity)
        new = self.evaluate(proposal,self.time+dt)
        proposal.velocity = half_velocity+.5*dt*new["force"]/self.node_mass[:,None]
        if not np.all(np.isfinite(proposal.velocity)):
            raise ValueError("Nonfinite velocity proposal; accepted state retained")
        self.state = proposal
        self.time += dt;self.steps += 1
        self.source_work += .5*dt*(old["source_power"]+new["source_power"])

    def advance(self, end_time, *, max_steps=100000, on_step=None):
        if not np.isfinite(end_time) or end_time < self.time:
            raise ValueError("Finite end time must not precede accepted time")
        count = 0
        while end_time-self.time > 8*np.finfo(float).eps*max(abs(end_time),1.):
            if count >= max_steps:
                raise RuntimeError("Step budget exceeded; accepted state retained")
            self.step(min(self.timestep_limit(),end_time-self.time));count += 1
            if on_step is not None:
                on_step(self)
        return self.state

    def diagnostics(self):
        e = self.evaluate();d = self.domain
        kinetic = float(.5*np.sum(self.node_mass[:,None]*self.state.velocity**2))
        energy = kinetic+e["internal_energy"]+e["magnetic"].energy
        boundary = self.state.coordinates[d.boundary]
        result = {"time_s":self.time,"steps":self.steps,"mass_kg":float(self.cell_mass.sum()),
                  "volume_m3":float(e["volume"].sum()),"kinetic_energy_j":kinetic,
                  "internal_energy_j":e["internal_energy"],"induced_field_energy_j":e["magnetic"].energy,
                  "response_energy_j":energy,"source_work_j":self.source_work,
                  "vacuum_relative_residual":e["magnetic"].vacuum_residual,
                  "magnetic_quadrature_difference_j":e["magnetic"].physical_energy-e["magnetic"].energy,
                  "minimum_density_kg_m3":float(np.min(self.cell_mass/e["volume"])),
                  "minimum_pressure_pa":float(e["pressure"].min()),
                  "maximum_speed_m_s":float(np.linalg.norm(self.state.velocity,axis=-1).max()),
                  "maximum_boundary_radius_m":float(boundary[:,0].max()),
                  "minimum_boundary_z_m":float(boundary[:,1].min()),"maximum_boundary_z_m":float(boundary[:,1].max()),
                  "axial_momentum_kg_m_s":float(self.node_mass@self.state.velocity[:,1]),
                  "minimum_jacobian_m2":min(float(b.jacobian.min()) for b in e["magnetic"].elements)}
        if hasattr(self,"initial"):
            result["response_energy_residual_j"] = energy-self.initial["response_energy_j"]-self.source_work
        return result
