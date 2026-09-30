"""Atomic density/two-entropy integration with explicitly prescribed bulk flow."""
from dataclasses import dataclass
import numpy as np
from .variable_compatible_current import VariableCompatibleCurrent
from .compatible_log_thermal import CompatibleLogThermal
from .compatible_particle_flux import CompatibleParticleFlux
from .compatible_bulk_thermal import CompatibleBulkThermal


@dataclass(frozen=True)
class BulkEntropyState:
    density: np.ndarray
    electron: np.ndarray
    ion: np.ndarray


@dataclass(frozen=True)
class BulkEntropyStep:
    state: BulkEntropyState
    thermal_work: np.ndarray
    max_scalar_residual: float


class PrescribedBulkTransport:
    """Verification driver: momentum supplied by a callback, never advanced here."""
    def __init__(self,geometry,ion_mass,momentum,*,gamma_e=5/3,gamma_i=1.4):
        self.geometry=geometry;self.mesh=geometry.mesh;self.ion_mass=ion_mass
        self.momentum=momentum;self.gammas=(gamma_e,gamma_i)

    def evaluate(self,state):
        model=VariableCompatibleCurrent(self.mesh,state.density,self.ion_mass,geometry=self.geometry)
        electron=CompatibleLogThermal(model,gamma=self.gammas[0])
        ion=CompatibleLogThermal(model,gamma=self.gammas[1],scalar=electron.scalar)
        fields=[t.fields(s) for t,s in zip((electron,ion),(state.electron,state.ion))]
        return model,(electron,ion),fields

    def derivative(self,state,time):
        model,species,_=self.evaluate(state);flux=CompatibleParticleFlux(model.density_map)
        faces=flux.from_momentum(self.momentum(model,time))
        rates=[CompatibleBulkThermal(t,flux).rate(s,faces)
               for t,s in zip(species,(state.electron,state.ion))]
        work=np.array([np.sum(faces.radial*r.face_work.radial)+np.sum(faces.axial*r.face_work.axial) for r in rates])
        residual=max(r.residuals[k] for r in rates for k in ('entropy_mass','thermal_projection'))
        return BulkEntropyState(rates[0].density,rates[0].entropy,rates[1].entropy),work,residual

    @staticmethod
    def combine(state,rate,step):
        return BulkEntropyState(*(np.asarray(a)+step*np.asarray(b) for a,b in
            zip((state.density,state.electron,state.ion),(rate.density,rate.electron,rate.ion))))

    def step(self,state,time,dt):
        if not np.isfinite([time,dt]).all() or dt<=0:raise ValueError('Finite time and positive step required')
        k1,w1,r1=self.derivative(state,time)
        k2,w2,r2=self.derivative(self.combine(state,k1,.5*dt),time+.5*dt)
        k3,w3,r3=self.derivative(self.combine(state,k2,.5*dt),time+.5*dt)
        k4,w4,r4=self.derivative(self.combine(state,k3,dt),time+dt)
        fields=[(getattr(k1,key)+2*getattr(k2,key)+2*getattr(k3,key)+getattr(k4,key))/6
                for key in ('density','electron','ion')]
        result=self.combine(state,BulkEntropyState(*fields),dt)
        self.evaluate(result)
        return BulkEntropyStep(result,dt*(w1+2*w2+2*w3+w4)/6,max(r1,r2,r3,r4))
