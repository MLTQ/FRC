"""Instantaneous finite-mass induction with moving-density and bulk-vorticity terms."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E
from .inertial_induction import InertialInduction,InductionRate
from .dual_kinematics import velocity
from .current_thermal_flux import current_particle_energy


@dataclass(frozen=True)
class MovingCurrentRate:
    induction: InductionRate
    electric_full: np.ndarray
    potential: np.ndarray


def nodal_drift(model,current):
    """Same explicit nodal current traces as the frozen fixed-n kinetic potential."""
    mesh=model.mesh;jr,jt,jz=model.geometry.split_electric(current)
    jr=jr.reshape(mesh.nr,mesh.nz);jz=jz.reshape(mesh.nr+1,mesh.nz)
    radial=np.zeros_like(jz);dx=mesh.radial_measure
    radial[1:-1]=mesh.rf[1:-1,None]*(dx[:-1,None]*jr[:-1]+dx[1:,None]*jr[1:])/(dx[:-1,None]+dx[1:,None])
    theta=np.pad(jt.reshape(mesh.nr-1,mesh.nz),((1,1),(0,0)))
    axial=.5*(jz+np.roll(jz,1,axis=1))
    return np.stack((radial,theta,axial),axis=-1)/(E*model.electron.density[...,None])


class MovingCurrent:
    """A moving-system rate component; density/velocity evolution is supplied by its caller."""
    def __init__(self,model,**solver_options):
        self.model=model;self.solver=InertialInduction(model,**solver_options)

    def terms(self,psi,g,ue,ui,regular_velocity,density_rate):
        model=self.model;geo=model.geometry;mesh=model.mesh
        if np.shape(psi)!=(mesh.nr+1,mesh.nz+1) or np.shape(g)!=(mesh.nr,mesh.nz) or not np.isfinite(psi).all() or not np.isfinite(g).all():
            raise ValueError('Finite corner flux and toroidal cell field required')
        scale=max(float(np.max(abs(psi))),np.finfo(float).tiny)
        if np.max(abs(psi[0]))>1e-13*scale or np.ptp(psi[-1])>1e-13*scale or np.max(abs(psi[:,0]-psi[:,-1]))>1e-13*scale:
            raise ValueError('Regular axis, fixed wall and periodic magnetic flux required')
        regular_velocity=np.asarray(regular_velocity,float)
        v,vorticity=velocity(geo,regular_velocity)
        if np.any(regular_velocity[-1,:,0]!=0):raise ValueError('Impermeable radial bulk boundary required')
        nt=model.electron.geometry.scalar(density_rate,'density time derivative')
        ntrate=geo.scalar_samples(nt)/model.density
        packed=geo.pack(psi,g);current=model.current(packed)
        w=geo.sample_electric(current)/(E*model.density[...,None]);star=model.star_field(packed,current)
        pressure,potential=model.pressure(ue,ui)
        terms=dict(hall=model.delta*np.cross(w,star),bulk=-np.cross(v,star),
                   vorticity=-model.reduced_mass/E*np.cross(w,vorticity),
                   density=-model.reduced_mass/E*w*ntrate[...,None])
        weak={k:model.mass@geo.project_electric(value) for k,value in terms.items()}
        weak['pressure']=model.mass@pressure
        vn=np.asarray(regular_velocity)*np.stack(np.broadcast_arrays(mesh.rf[:,None],mesh.rf[:,None],np.ones((mesh.nr+1,mesh.nz))),axis=-1)
        potential=potential+model.reduced_mass/E*np.sum(vn*nodal_drift(model,current),axis=-1)
        potential-=model.delta/E*current_particle_energy(model,current)
        return weak,potential

    def rate(self,psi,g,ue,ui,regular_velocity,density_rate,*,electric_source=None):
        weak,potential=self.terms(psi,g,ue,ui,regular_velocity,density_rate)
        forcing=sum(weak.values())
        if electric_source is not None:
            source=np.asarray(electric_source,float)
            if source.shape!=self.model.geometry.r.shape+(3,) or not np.isfinite(source).all():
                raise ValueError('Finite externally specified electric source at physical quadrature required')
            forcing+=self.model.mass@self.model.geometry.project_electric(source)
        solved=self.solver.solve(forcing)
        return MovingCurrentRate(solved,solved.electric+self.model.gradient(potential),potential)
