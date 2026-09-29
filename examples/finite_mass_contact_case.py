"""Projection and independent contact errors for the finite-mass state."""
import numpy as np
from plasma_bench.finite_mass_state import FiniteMassParameters,from_primitive,evaluate
from plasma_bench.dual_mass import DualMass


def averages(mesh,ref,time=0.,*,rate=False,order=12):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    fields=ref.fields(mesh.z[:,None]+.5*mesh.dz*nodes,time,rate=rate)
    return {key:np.broadcast_to(.5*np.sum(value*weights,axis=1)[None],(mesh.nr,mesh.nz)).copy() for key,value in fields.items()}


def initialize(mesh,ref):
    par=FiniteMassParameters(ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    exact=averages(mesh,ref);nodes=ref.fields(mesh.zf[:-1]);pe=np.broadcast_to(nodes['electron_pressure'][None],(mesh.nr+1,mesh.nz)).copy()
    # Complement the actual Q1 electron pressure average so the discrete initial
    # total pressure is constant. The projection's pi/energy error is recorded.
    pi=ref.total_pressure-DualMass(mesh).primal_average(pe)
    psi=.5*ref.bias*mesh.corners[0]**2
    return from_primitive(mesh,par,psi,0.,exact['density'],pe,pi,vz=ref.velocity),par


def norms(mesh,ref,state,par,time,*,order=12):
    ev=evaluate(mesh,state,par);exact=averages(mesh,ref,time,order=order)
    pe=ev.model.electron.pressure(state.electron_entropy);pi=(par.gamma_i-1)*ev.components['ion']
    nodal=ref.fields(mesh.zf[:-1],time)
    items=dict(density=(state.conserved[...,0],exact['density'],mesh.volumes),
        momentum=(state.conserved[...,2],exact['momentum'],mesh.volumes),
        total_energy=(state.conserved[...,3],exact['energy'],mesh.volumes),
        electron_pressure=(pe,np.broadcast_to(nodal['electron_pressure'][None],pe.shape),ev.model.mapping.volumes),
        ion_pressure=(pi,exact['ion_pressure'],mesh.volumes),
        total_pressure=(ev.primitive[...,3],np.full_like(pi,ref.total_pressure),mesh.volumes),
        velocity=(ev.primitive[...,2],np.full_like(pi,ref.velocity),mesh.volumes))
    return {key:dict(relative_L2=float(np.sqrt(np.sum(v*(a-b)**2)/np.sum(v*b*b))),
        rms_error=float(np.sqrt(np.sum(v*(a-b)**2)/np.sum(v))),maximum_error=float(abs(a-b).max())) for key,(a,b,v) in items.items()}
