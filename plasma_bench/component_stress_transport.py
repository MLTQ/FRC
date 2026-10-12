"""Evaluate original drift-stress entries only as the adjoint requests them."""
import numpy as np
from scipy.constants import elementary_charge as E
from .scheduled_current_transport import ScheduledCurrentTransport


class ComponentDriftStress:
    """Small adapter for the existing stress_rate component-index contract."""
    def __init__(self, factor, current, density):
        self.factor = factor; self.current = current; self.density = density

    def __getitem__(self, key):
        if (not isinstance(key, tuple) or len(key) != 3 or key[0] is not Ellipsis
                or any(type(k) is not int or not 0 <= k < 3 for k in key[1:])):
            raise IndexError('Only explicit physical stress components are supported')
        i, j = key[1:]
        return self.factor*self.current[..., i]*self.current[..., j]/self.density


class ComponentStressTransport(ScheduledCurrentTransport):
    def pair(self,state,current,density_rate):
        m=self.model;g=m.geometry
        if not np.array_equal(state[...,0],m.density_map.point_values):raise ValueError('Same actual fluid/current density required')
        velocity=m.mapping.velocity(state);j=self.current.sample(current);n=m.density;w=j/(E*n[...,None])
        stress=ComponentDriftStress(m.reduced_mass/E**2,j,n)
        bulk=self.fluid.stress_rate(state,stress)
        del stress
        grad=self.fluid.velocity_gradient(state)
        stretching=m.reduced_mass/E**2*self.current.adjoint(np.einsum('...ij,...j->...i',grad,j)/n[...,None])
        del grad, j
        face_v=[];face_w=[]
        for radial,surface,maps,fluids in self.faces:
            nl,nr=(fm.scalar(state[...,0]) for fm in fluids);vl,vr=(fm.velocity(state) for fm in fluids);jl,jr=(cm.sample(current) for cm in maps)
            face_v.append((vl,vr,nl,nr));face_w.append((jl/(E*nl[...,None]),jr/(E*nr[...,None]),nl,nr))
            if radial:
                traction=m.reduced_mass/(2*E**2)*(jl*(jl[...,0]/nl)[...,None]+jr*(jr[...,0]/nr)[...,None])
                bulk-=fluids[0].force_rate(state,traction);bulk+=fluids[1].force_rate(state,traction)
                jump=vl-vr
                stretching-=m.reduced_mass/(2*E**2)*(maps[0].adjoint(jump*(jl[...,0]/nl)[...,None])+maps[1].adjoint(jump*(jr[...,0]/nr)[...,None]))
                del traction, jump
        advection=self.advection(current,velocity,density_rate,face_v)
        del velocity, face_v
        self_drift=-m.delta*self.advection(current,w,np.zeros_like(density_rate),face_w)
        del w, face_w
        result=dict(fluid_rate=bulk,advection=advection,stretching=stretching,self_drift=self_drift,forcing=advection+stretching+self_drift)
        if any(not np.isfinite(v).all() for v in result.values()):raise ValueError('Nonfinite physical transport')
        return result
