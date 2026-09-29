"""Unsplit continuity, momentum, species entropy and finite-inertia induction."""
from dataclasses import dataclass
import numpy as np
from .finite_mass_state import FiniteMassState,evaluate
from .finite_mass_flux import scalar_faces,rusanov
from .finite_mass_transport import add_transport
from .primal_velocity import regular_velocity
from .toroidal_reconstruction import ToroidalReconstruction
from .toroidal_state import radius_squared
from .bulk_entropy import entropy_rate
from .moving_current import MovingCurrent


@dataclass(frozen=True)
class FiniteMassSources:
    fluid: np.ndarray
    electron_entropy: np.ndarray
    electric: np.ndarray


@dataclass
class FiniteMassRate:
    state: FiniteMassState
    electric: np.ndarray
    radial_flux: np.ndarray
    axial_flux: np.ndarray
    fluid_source: np.ndarray
    evaluation: object
    number_density_rate: np.ndarray
    regular_velocity: np.ndarray
    linear_residual: float
    linear_iterations: int
    linear_refinements: int


def rhs(mesh,state,parameters,*,time=0.,source=None):
    ev=evaluate(mesh,state,parameters);model=ev.model;mapping=model.mapping;q=ev.primitive
    reconstruct=ToroidalReconstruction(mesh,q,order=parameters.spatial_order)
    left,right,down,up=reconstruct.faces()
    fraction=(parameters.gamma_e-1)*ev.components['electron']/q[...,3]
    fractions=scalar_faces(mesh,fraction,parameters.spatial_order,upper=1.)
    kinetic=scalar_faces(mesh,ev.components['current'],parameters.spatial_order)
    br,bz=mesh.face_fields(state.psi)
    fr=rusanov(left,right,*fractions[:2],*kinetic[:2],br,mesh.rf[:,None],0,parameters)
    fz=rusanov(down,up,*fractions[2:],*kinetic[2:],bz,np.sqrt(radius_squared(mesh)),1,parameters)
    fr*=2*np.pi*mesh.rf[:,None,None]*mesh.dz
    fz*=2*np.pi*mesh.radial_measure[:,None,None]
    mass_rate=-(np.diff(fr[...,0],axis=0)+np.diff(fz[...,0],axis=1))/mesh.volumes
    nt=mapping.density(mass_rate)/parameters.masses.total
    particle=mapping.fluxes(fr[...,0]/parameters.masses.total,fz[...,0]/parameters.masses.total)
    regular=regular_velocity(mapping,state.conserved)
    forcing=source(time,ev) if source is not None else None
    fluid_source=np.zeros_like(state.conserved);electron_source=np.zeros_like(state.electron_entropy)
    if forcing is not None:
        shapes=(fluid_source.shape,electron_source.shape,model.geometry.r.shape+(3,))
        if any(np.shape(a)!=shape or not np.isfinite(a).all() for a,shape in zip((forcing.fluid,forcing.electron_entropy,forcing.electric),shapes)):
            raise ValueError('Finite explicit momentum/energy, entropy and electric sources required')
        if np.any(forcing.fluid[...,0]!=0):raise ValueError('Particle creation/ionization source is not implemented')
        fluid_source=forcing.fluid.copy();electron_source=forcing.electron_entropy
    induction=MovingCurrent(model,tolerance=parameters.linear_tolerance,max_iterations=parameters.max_iterations).rate(
        state.psi,state.g,state.electron_entropy,ev.ion_entropy,regular,nt,
        electric_source=None if forcing is None else forcing.electric)
    hoop=add_transport(state,ev,regular,induction.electric_full,fr,fz)
    du=-(np.diff(fr,axis=0)+np.diff(fz,axis=1))/mesh.volumes[...,None]+fluid_source
    du[...,1]+=reconstruct.radial_source()+hoop
    electron=entropy_rate(mapping,model.electron.density,state.electron_entropy,*particle,order=parameters.spatial_order)
    electron+=model.beta*model.electron.rhs(state.g,state.electron_entropy)[1]+electron_source
    p,g=model.geometry.unpack(induction.induction.magnetic)
    derivative=FiniteMassState(du,p,g,electron)
    return FiniteMassRate(derivative,induction.electric_full,fr,fz,fluid_source,ev,nt,regular,
        induction.induction.residual,induction.induction.iterations,induction.induction.refinements)


def pressure_rates(state,rate,parameters):
    """Differentiate the stated energy partition, including physical J²/n and n_t."""
    ev=rate.evaluation;model=ev.model;geo=model.geometry;mesh=model.mesh;q=ev.primitive
    du=rate.state.conserved;vr,vz,omega=q[...,1],q[...,2],q[...,6]
    bulk=vr*du[...,1]+vz*du[...,2]+omega*du[...,4]
    bulk-=.5*(vr*vr+vz*vz+radius_squared(mesh)*omega*omega)*du[...,0]
    b=geo.magnetic_samples(state.psi,state.g);bd=geo.magnetic_samples(rate.state.psi,rate.state.g)
    from scipy.constants import mu_0 as MU0
    magnetic=np.sum(geo.weight*np.sum(b*bd,axis=-1),axis=(2,3))/(MU0*mesh.volumes)
    j=geo.sample_electric(ev.current);jd=geo.sample_electric(model.current(geo.pack(rate.state.psi,rate.state.g)))
    nt=geo.scalar_samples(rate.number_density_rate)
    current=np.sum(geo.weight*model.inertia_coefficient*(np.sum(j*jd,axis=-1)-.5*np.sum(j*j,axis=-1)*nt/model.density),axis=(2,3))/mesh.volumes
    pe=model.electron.pressure(state.electron_entropy)
    electron=model.electron.a*rate.state.electron_entropy+(parameters.gamma_e-1)*pe*rate.number_density_rate/model.electron.density
    ion=(parameters.gamma_i-1)*(du[...,3]-bulk-magnetic-current-model.mapping.primal_average(electron)/(parameters.gamma_e-1))
    return electron,ion
