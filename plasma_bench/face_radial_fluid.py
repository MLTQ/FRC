"""Self-consistent fluid with one common radial transport/kinetic momentum field."""
from dataclasses import dataclass,fields
import numpy as np
from scipy.constants import electron_mass
from .variable_compatible_current import VariableCompatibleCurrent
from .compatible_log_thermal import CompatibleLogThermal
from .compatible_particle_flux import CompatibleParticleFlux,ParticleFaces
from .compatible_bulk_thermal import CompatibleBulkThermal
from .face_radial_kinematics import FaceRadialKinematics
from .face_radial_flux import FaceRadialFlux


@dataclass(frozen=True)
class FaceRadialState:
    density: np.ndarray
    radial: np.ndarray
    momentum: np.ndarray
    electron: np.ndarray
    ion: np.ndarray
    energy_account: np.ndarray


class FaceRadialFluid:
    def __init__(self,geometry,ion_mass,*,gamma_e=5/3,gamma_i=1.4):
        self.geometry=geometry;self.mesh=geometry.mesh;self.ion_mass=ion_mass;self.mass=ion_mass+electron_mass
        self.gammas=(gamma_e,gamma_i);self.kinematics=FaceRadialKinematics(geometry,self.mass)
        self.transport=FaceRadialFlux(geometry,self.mass,self.kinematics)

    def evaluate(self,state):
        if np.shape(state.energy_account)!=(self.mesh.nr,self.mesh.nz) or not np.isfinite(state.energy_account).all():
            raise ValueError('Finite conservative diagnostic energy account required')
        model=VariableCompatibleCurrent(self.mesh,state.density,self.ion_mass,geometry=self.geometry)
        e=CompatibleLogThermal(model,gamma=self.gammas[0]);i=CompatibleLogThermal(model,gamma=self.gammas[1],scalar=e.scalar)
        f=[t.fields(s) for t,s in zip((e,i),(state.electron,state.ion))]
        kinetic=self.kinematics.fields(model.density_map,state.radial,state.momentum)
        energy=model.average(kinetic['kinetic']+f[0]['internal']+f[1]['internal'])
        return model,(e,i),f,kinetic,energy

    def initialize(self,density,radial,momentum,electron,ion):
        s=FaceRadialState(*(np.asarray(v,float).copy() for v in (density,radial,momentum,electron,ion)),np.zeros_like(density))
        return FaceRadialState(s.density,s.radial,s.momentum,s.electron,s.ion,self.evaluate(s)[-1])

    def particle_faces(self,flux,state):
        return ParticleFaces(self.kinematics.full_faces(state.radial),np.asarray(state.momentum[...,0]@flux.axial_face.T))

    def radial_kinetic_work(self,model,flux,state,mechanical,kinetic,gradient):
        g=self.geometry;m=self.mesh;density=model.density_map;v=kinetic['velocity']
        cell=self.mass*self.kinematics.cell_adjoint(density,v)/m.volumes[...,None]
        jumps=np.zeros((m.nr+1,m.nz,2));jumps[1:-1]=np.diff(cell,axis=0)
        point_load=np.sum(np.repeat(jumps,g.order,axis=1)*mechanical['radial_specific_momentum'],axis=-1)
        weights=np.tile(np.pi*m.dz*self.transport.gauss_weights,m.nz)
        cell_advection=np.asarray((point_load*weights[None,:])@density.axial)
        # Radial projected convection contributes -M vr d_r(vr) Gamma_r
        # and M vr^2 n_t. Combine its density work with true K_n.
        density_load=density.adjoint(.5*self.mass*(v[...,0]**2-v[...,1]**2-v[...,2]**2))
        vector=np.zeros_like(v);vector[...,0]=-self.mass*v[...,0]*gradient['radial'][...,0]
        return (cell_advection+flux.density_adjoint(density_load).radial+flux.vector_adjoint(vector).radial)

    def derivative(self,state,time=0.):
        model,species,_,kinetic,_=self.evaluate(state);flux=CompatibleParticleFlux(model.density_map)
        faces=self.particle_faces(flux,state);states=(state.electron,state.ion)
        thermal=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,states)]
        mechanical=self.transport.rate(model,state.radial,state.momentum,species,states)
        nt=flux.density_rate(faces);particle=flux.sample(faces)
        gradient=self.kinematics.gradients(model.density_map,state.radial,state.momentum);v=kinetic['velocity']
        radial_advection=(-particle[...,0]*gradient['radial'][...,0]-particle[...,2]*gradient['axial'][...,0]
                          +v[...,0]*model.density_map.sample(nt))
        # Keep the physical centrifugal force. Cancelling the discrete cell-
        # momentum convective energy defect here would change the Euler force.
        centrifugal=model.density*v[...,1]**2/self.geometry.r
        work=sum(t.face_work.radial for t in thermal)
        qt,residual=self.kinematics.solve_rate(model.density_map,radial_advection+centrifugal,work)
        mt=mechanical['momentum'];rate=FaceRadialState(nt,qt,mt,thermal[0].entropy,thermal[1].entropy,mechanical['energy_account'])
        kt=self.kinematics.energy_rate(model.density_map,state.radial,state.momentum,nt,qt,mt)['total']
        ht=[t.energy_rate(s,k.entropy,nt)['total'] for t,s,k in zip(species,states,thermal)]
        expected=flux.integrated(faces)
        continuity=max(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300)
                       for a,b in zip(expected,(mechanical['radial_particles'],mechanical['axial_particles'])))
        return rate,dict(local_energy_defect=kt+sum(ht)-mechanical['energy_account'],kinetic_rate=kt,thermal_rates=ht,
            radial_metric_residual=residual,continuity_face_agreement=continuity,
            scalar_residual=max([residual]+[k.residuals[key] for k in thermal for key in ('entropy_mass','thermal_projection')]))

    @staticmethod
    def combine(state,rate,step):
        return FaceRadialState(*(getattr(state,f.name)+step*getattr(rate,f.name) for f in fields(state)))

    def step(self,state,time,dt):
        if not np.isfinite([time,dt]).all() or dt<=0:raise ValueError('Finite time and positive step required')
        k1,d1=self.derivative(state,time);k2,d2=self.derivative(self.combine(state,k1,.5*dt),time+.5*dt)
        k3,d3=self.derivative(self.combine(state,k2,.5*dt),time+.5*dt);k4,d4=self.derivative(self.combine(state,k3,dt),time+dt)
        rate=FaceRadialState(*((getattr(k1,f.name)+2*getattr(k2,f.name)+2*getattr(k3,f.name)+getattr(k4,f.name))/6 for f in fields(state)))
        result=self.combine(state,rate,dt);self.evaluate(result)
        ledger=dt*(d1['local_energy_defect']+2*d2['local_energy_defect']+2*d3['local_energy_defect']+d4['local_energy_defect'])/6
        return result,dict(integrated_spatial_energy_defect=ledger,
            scalar_residual=max(d['scalar_residual'] for d in (d1,d2,d3,d4)),
            continuity_face_agreement=max(d['continuity_face_agreement'] for d in (d1,d2,d3,d4)))
