"""Self-consistent bulk motion and two physical log entropies in a closed cylinder."""
from dataclasses import dataclass,fields
import numpy as np
from scipy.constants import electron_mass
from .variable_compatible_current import VariableCompatibleCurrent
from .compatible_log_thermal import CompatibleLogThermal
from .compatible_particle_flux import CompatibleParticleFlux
from .compatible_bulk_thermal import CompatibleBulkThermal
from .conservative_kinematics import ConservativeKinematics
from .log_fluid_flux import LogFluidFlux


@dataclass(frozen=True)
class LogFluidState:
    density: np.ndarray
    momentum: np.ndarray
    electron: np.ndarray
    ion: np.ndarray
    energy_account: np.ndarray


class LogFluid:
    """Smooth hydrodynamic component; magnetic/current dynamics are not advanced."""
    def __init__(self,geometry,ion_mass,*,gamma_e=5/3,gamma_i=1.4):
        self.geometry=geometry;self.mesh=geometry.mesh;self.ion_mass=ion_mass
        self.gammas=(gamma_e,gamma_i);self.mass=ion_mass+electron_mass
        self.kinematics=ConservativeKinematics(geometry,self.mass)
        self.transport=LogFluidFlux(geometry,self.mass)

    def evaluate(self,state):
        if np.shape(state.energy_account)!=(self.mesh.nr,self.mesh.nz) or not np.isfinite(state.energy_account).all():
            raise ValueError('Finite conservative diagnostic energy account required')
        model=VariableCompatibleCurrent(self.mesh,state.density,self.ion_mass,geometry=self.geometry)
        electron=CompatibleLogThermal(model,gamma=self.gammas[0])
        ion=CompatibleLogThermal(model,gamma=self.gammas[1],scalar=electron.scalar)
        f=[t.fields(s) for t,s in zip((electron,ion),(state.electron,state.ion))]
        kinetic=self.kinematics.fields(model.density_map,state.momentum)
        energy=model.average(kinetic['kinetic']+f[0]['internal']+f[1]['internal'])
        return model,(electron,ion),f,kinetic,energy

    def initialize(self,density,momentum,electron,ion):
        state=LogFluidState(*(np.asarray(v,float).copy() for v in (density,momentum,electron,ion)),np.zeros_like(density))
        energy=self.evaluate(state)[-1]
        return LogFluidState(state.density,state.momentum,state.electron,state.ion,energy)

    def derivative(self,state,time=0.):
        model,species,_,kinetic,_=self.evaluate(state);flux=CompatibleParticleFlux(model.density_map)
        faces=flux.from_momentum(state.momentum[...,:2]);states=(state.electron,state.ion)
        thermal=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,states)]
        mechanical=self.transport.rate(model,self.kinematics,state.momentum,species,states)
        nt=flux.density_rate(faces);mt=mechanical['momentum']
        rate=LogFluidState(nt,mt,thermal[0].entropy,thermal[1].entropy,mechanical['energy_account'])
        kt=self.kinematics.energy_rate(model.density_map,state.momentum,nt,mt)['total']
        ht=[t.energy_rate(s,k.entropy,nt)['total'] for t,s,k in zip(species,states,thermal)]
        local_defect=kt+sum(ht)-mechanical['energy_account']
        expected=flux.integrated(faces)
        continuity_error=max(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300)
            for a,b in zip(expected,(mechanical['radial_particles'],mechanical['axial_particles'])))
        diagnostics=dict(local_energy_defect=local_defect,kinetic_rate=kt,thermal_rates=ht,
            continuity_face_agreement=continuity_error,
            scalar_residual=max(k.residuals[key] for k in thermal for key in ('entropy_mass','thermal_projection')))
        return rate,diagnostics

    @staticmethod
    def combine(state,rate,step):
        return LogFluidState(*(getattr(state,f.name)+step*getattr(rate,f.name) for f in fields(state)))

    def step(self,state,time,dt):
        if not np.isfinite([time,dt]).all() or dt<=0:raise ValueError('Finite time and positive step required')
        k1,d1=self.derivative(state,time)
        k2,d2=self.derivative(self.combine(state,k1,.5*dt),time+.5*dt)
        k3,d3=self.derivative(self.combine(state,k2,.5*dt),time+.5*dt)
        k4,d4=self.derivative(self.combine(state,k3,dt),time+dt)
        rate=LogFluidState(*((getattr(k1,f.name)+2*getattr(k2,f.name)+2*getattr(k3,f.name)+getattr(k4,f.name))/6 for f in fields(state)))
        result=self.combine(state,rate,dt);self.evaluate(result)
        ledger=dt*(d1['local_energy_defect']+2*d2['local_energy_defect']+2*d3['local_energy_defect']+d4['local_energy_defect'])/6
        return result,dict(integrated_spatial_energy_defect=ledger,
            scalar_residual=max(d['scalar_residual'] for d in (d1,d2,d3,d4)),
            continuity_face_agreement=max(d['continuity_face_agreement'] for d in (d1,d2,d3,d4)))
