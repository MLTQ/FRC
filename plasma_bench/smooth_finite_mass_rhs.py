"""Unsplit finite-mass derivative with shared smooth higher-order bulk fluxes."""
from dataclasses import dataclass
import numpy as np
from .finite_mass_state import FiniteMassState,evaluate
from .finite_mass_rhs import FiniteMassRate
from .smooth_finite_mass_flux import bulk_fluxes
from .finite_mass_transport import add_transport
from .primal_velocity import regular_velocity
from .bulk_entropy import entropy_rate as limited_entropy_rate
from .smooth_bulk_entropy import entropy_rate
from .kernel_induction import KernelMovingCurrent


@dataclass
class SmoothFiniteMassRate(FiniteMassRate):
    induction_components:dict


def rhs(mesh,state,parameters,*,time=0.,source=None,quadrature=5,entropy_transport='smooth'):
    if entropy_transport not in ('smooth','limited'):raise ValueError('Explicit smooth or limited entropy transport required')
    ev=evaluate(mesh,state,parameters);model=ev.model;mapping=model.mapping
    fr,fz,bulk_hoop=bulk_fluxes(mesh,state,ev,parameters,quadrature=quadrature)
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
    induction=KernelMovingCurrent(model,tolerance=parameters.linear_tolerance,max_iterations=parameters.max_iterations).rate(
        state.psi,state.g,state.electron_entropy,ev.ion_entropy,regular,nt,
        electric_source=None if forcing is None else forcing.electric)
    hoop=add_transport(state,ev,regular,induction.electric_full,fr,fz)
    du=-(np.diff(fr,axis=0)+np.diff(fz,axis=1))/mesh.volumes[...,None]+fluid_source
    du[...,1]+=bulk_hoop+hoop
    if entropy_transport=='smooth':electron=entropy_rate(mapping,model.electron.density,state.electron_entropy,*particle)
    else:electron=limited_entropy_rate(mapping,model.electron.density,state.electron_entropy,*particle,order=parameters.spatial_order)
    electron+=model.beta*model.electron.rhs(state.g,state.electron_entropy)[1]+electron_source
    p,g=model.geometry.unpack(induction.induction.magnetic)
    return SmoothFiniteMassRate(FiniteMassState(du,p,g,electron),induction.electric_full,fr,fz,fluid_source,ev,nt,regular,
        induction.induction.residual,induction.induction.iterations,induction.induction.refinements,induction.induction.residual_components)
