"""Experimental fluid variant with physical radial thermal/kinetic work pairing."""
import numpy as np
from .log_fluid import LogFluid,LogFluidState
from .compatible_particle_flux import CompatibleParticleFlux
from .compatible_bulk_thermal import CompatibleBulkThermal
from .radial_pressure_metric import RadialPressureMetric


class RadialPressureFluid(LogFluid):
    """Diagnostic radial work pair; retains conservative axial/angular transport."""

    def radial_kinetic_work(self,model,flux,state):
        g=self.geometry;m=self.mesh;density=model.density_map
        gradients=self.kinematics.energy_gradients(density,state.momentum)
        cell=gradients['momentum']/m.volumes[...,None]
        jumps=np.zeros((m.nr+1,m.nz,3));jumps[1:-1]=np.diff(cell,axis=0)
        maps=self.transport.rfaces
        n=np.asarray(maps[0]@state.density@density.axial.T)
        momenta=np.stack([np.asarray(maps[p]@state.momentum[...,k]@density.axial.T)
                          for p,k in ((1,0),(0,1),(2,2))],axis=-1)
        momenta[[0,-1],...,0]=0.
        point_load=np.sum(np.repeat(jumps,g.order,axis=1)*momenta/n[...,None],axis=-1)
        weights=np.tile(np.pi*m.dz*self.transport.gauss_weights,m.nz)
        advection=np.asarray((point_load*weights[None,:])@density.axial)
        return advection+flux.density_adjoint(gradients['density']).radial

    def derivative(self,state,time=0.):
        model,species,_,_,_=self.evaluate(state);flux=CompatibleParticleFlux(model.density_map)
        faces=flux.from_momentum(state.momentum[...,:2]);states=(state.electron,state.ion)
        thermal=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,states)]
        mechanical=self.transport.rate(model,self.kinematics,state.momentum,species,states)
        pressure_only=self.transport.rate(model,self.kinematics,np.zeros_like(state.momentum),species,states)
        if not hasattr(self,'radial_pressure'):self.radial_pressure=RadialPressureMetric(self.kinematics,model.density_map)
        work=sum(t.face_work.radial for t in thermal)+self.radial_kinetic_work(model,flux,state)
        radial,residual=self.radial_pressure.solve(model,flux,work)
        nt=flux.density_rate(faces);mt=mechanical['momentum'].copy()
        velocity=self.kinematics.fields(model.density_map,state.momentum)['velocity']
        centrifugal=model.average(model.density*velocity[...,1]**2/self.geometry.r)
        mt[...,0]+=radial-pressure_only['momentum'][...,0]-centrifugal
        rate=LogFluidState(nt,mt,thermal[0].entropy,thermal[1].entropy,mechanical['energy_account'])
        kt=self.kinematics.energy_rate(model.density_map,state.momentum,nt,mt)['total']
        ht=[t.energy_rate(s,k.entropy,nt)['total'] for t,s,k in zip(species,states,thermal)]
        expected=flux.integrated(faces)
        continuity_error=max(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300)
            for a,b in zip(expected,(mechanical['radial_particles'],mechanical['axial_particles'])))
        diagnostics=dict(local_energy_defect=kt+sum(ht)-mechanical['energy_account'],kinetic_rate=kt,thermal_rates=ht,
            continuity_face_agreement=continuity_error,radial_pressure_residual=residual,
            scalar_residual=max([residual]+[k.residuals[key] for k in thermal for key in ('entropy_mass','thermal_projection')]))
        return rate,diagnostics
