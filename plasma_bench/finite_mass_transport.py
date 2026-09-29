"""Full current energy and drift-stress fluxes through shared primal faces."""
import numpy as np
from scipy.constants import elementary_charge as E
from .compatible_poynting import energy_flux
from .dual_enthalpy import scalar_current_flux
from .moving_current import nodal_drift


def add_transport(state,evaluation,regular,electric,radial,axial):
    """Add integrated drift stress/energy and full Poynting power; return hoop source."""
    model=evaluation.model;mesh=model.mesh;geo=model.geometry;dual=model.electron.geometry
    w=nodal_drift(model,evaluation.current)
    v=regular*np.stack(np.broadcast_arrays(mesh.rf[:,None],mesh.rf[:,None],np.ones(dual.volumes.shape)),axis=-1)
    charge=dual.current_flux(state.g,0.)
    for column,value in ((1,w[...,0]),(2,w[...,2]),(4,mesh.rf[:,None]*w[...,1])):
        dr,dz=scalar_current_flux(dual,*charge,model.reduced_mass/E*value)
        radial[...,column]+=dr;axial[...,column]+=dz
    pe=model.electron.pressure(state.electron_entropy);pi=model.ion.pressure(evaluation.ion_entropy)
    he=model.electron.gamma/(model.electron.gamma-1)*pe
    hi=model.ion.gamma/(model.ion.gamma-1)*pi
    kinetic=.5*model.reduced_mass*np.sum(w*w,axis=-1)
    scalar=(model.reduced_mass*np.sum(v*w,axis=-1)-model.delta*kinetic
            +(model.alpha*hi-model.beta*he)/model.electron.density)/E
    dr,dz=scalar_current_flux(dual,*charge,scalar)
    sr,sz=energy_flux(geo,state.psi,state.g,electric)
    radial[...,3]+=dr+2*np.pi*mesh.rf[:,None]*mesh.dz*sr
    axial[...,3]+=dz+2*np.pi*mesh.radial_measure[:,None]*sz
    current=geo.sample_electric(evaluation.current)
    return np.sum(geo.weight*model.inertia_coefficient*current[...,1]**2/geo.r,axis=(2,3))/mesh.volumes
