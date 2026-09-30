"""Candidate shared current-advection/stress work with explicit face traces."""
import numpy as np
from scipy.constants import elementary_charge as E
from .physical_derivatives import CurrentDifferential,RegularDifferential,TensorSurface,fluid_map


class CurrentTransport:
    def __init__(self,model):
        self.model=model;g=model.geometry;self.current=CurrentDifferential(g);self.fluid=RegularDifferential(model.mapping);self.faces=[]
        for radial in (True,False):
            surface=TensorSurface(g,radial);count=g.mesh.nr if radial else g.mesh.nz
            right=np.arange(1,count) if radial else np.arange(count);left=right-1
            choices=({'radial_cells':left},{'radial_cells':right}) if radial else ({'axial_cells':left},{'axial_cells':right})
            currents=[CurrentDifferential(surface,**choice) for choice in choices]
            fluids=[fluid_map(model.mapping.fluid,surface,side if radial else None) for side in ('left','right')]
            self.faces.append((radial,surface,currents,fluids))

    def advection(self,current,velocity,density_rate,face_values):
        m=self.model;g=m.geometry;j=self.current.sample(current);gradient=self.current.gradient(current)
        adv=np.einsum('...ij,...j->...i',gradient,velocity);factor=m.reduced_mass/(2*E**2)
        result=factor*(self.current.adjoint(adv/m.density[...,None])-self.current.advective_adjoint(j/m.density[...,None],velocity))
        result+=factor*self.current.adjoint(j*(m.density_map.sample(density_rate)/m.density**2)[...,None])
        for (radial,surface,maps,fluids),(vl,vr,nl,nr) in zip(self.faces,face_values):
            normal=0 if radial else 2;jl,jr=(cm.sample(current) for cm in maps)
            speed=.5*(vl[...,normal]/nl+vr[...,normal]/nr)
            result+=factor*(maps[0].adjoint(speed[...,None]*jr)-maps[1].adjoint(speed[...,None]*jl))
        return result

    def pair(self,state,current,density_rate):
        m=self.model;g=m.geometry
        if not np.array_equal(state[...,0],m.density_map.point_values):raise ValueError('Same actual fluid/current density required')
        velocity=m.mapping.velocity(state);j=self.current.sample(current);n=m.density;w=j/(E*n[...,None]);grad=self.fluid.velocity_gradient(state)
        stress=m.reduced_mass/E**2*j[..., :,None]*j[...,None,:]/n[...,None,None]
        bulk=self.fluid.stress_rate(state,stress)
        stretching=m.reduced_mass/E**2*self.current.adjoint(np.einsum('...ij,...j->...i',grad,j)/n[...,None])
        face_v=[];face_w=[]
        for radial,surface,maps,fluids in self.faces:
            nl,nr=(fm.scalar(state[...,0]) for fm in fluids);vl,vr=(fm.velocity(state) for fm in fluids);jl,jr=(cm.sample(current) for cm in maps)
            face_v.append((vl,vr,nl,nr));face_w.append((jl/(E*nl[...,None]),jr/(E*nr[...,None]),nl,nr))
            if radial:
                traction=m.reduced_mass/(2*E**2)*(jl*(jl[...,0]/nl)[...,None]+jr*(jr[...,0]/nr)[...,None])
                bulk-=fluids[0].force_rate(state,traction);bulk+=fluids[1].force_rate(state,traction)
                jump=vl-vr
                stretching-=m.reduced_mass/(2*E**2)*(maps[0].adjoint(jump*(jl[...,0]/nl)[...,None])+maps[1].adjoint(jump*(jr[...,0]/nr)[...,None]))
        advection=self.advection(current,velocity,density_rate,face_v)
        self_drift=-m.delta*self.advection(current,w,np.zeros_like(density_rate),face_w)
        result=dict(fluid_rate=bulk,advection=advection,stretching=stretching,self_drift=self_drift,forcing=advection+stretching+self_drift)
        if any(not np.isfinite(v).all() for v in result.values()):raise ValueError('Nonfinite physical transport')
        return result
