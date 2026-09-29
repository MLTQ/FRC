"""Independent conservative-state and species errors of manufactured trajectories."""
import numpy as np
from plasma_bench.finite_mass_state import evaluate
from .finite_mass_manufactured import cell_average,point_conserved


def measure(mesh,ref,initial,state,parameters,start,end,*,order=8):
    exact0=cell_average(mesh,lambda r,z:point_conserved(ref,r,z,start),order)
    exact1=cell_average(mesh,lambda r,z:point_conserved(ref,r,z,end),order)
    ev=evaluate(mesh,state,parameters);ev0=evaluate(mesh,initial,parameters)
    items={key:(initial.conserved[...,i],state.conserved[...,i],exact0[...,i],exact1[...,i],mesh.volumes)
           for i,key in enumerate(('mass','radial_momentum','axial_momentum','total_energy','angular_momentum'))}
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij')
    items['electron_pressure']=(ev0.model.electron.pressure(initial.electron_entropy),ev.model.electron.pressure(state.electron_entropy),
        ref.fields(r,z,start)['pe'],ref.fields(r,z,end)['pe'],ev.model.mapping.volumes)
    items['ion_pressure']=((parameters.gamma_i-1)*ev0.components['ion'],(parameters.gamma_i-1)*ev.components['ion'],
        cell_average(mesh,lambda r,z:ref.fields(r,z,start)['pi'],order),
        cell_average(mesh,lambda r,z:ref.fields(r,z,end)['pi'],order),mesh.volumes)
    result={}
    for key,(a,b,x,y,v) in items.items():
        norm=lambda z:float(np.sqrt(np.sum(v*z*z)))
        error=b-y;change_error=(b-a)-(y-x);den=max(norm(y-x),np.finfo(float).tiny)
        result[key]=dict(relative_L2=norm(error)/norm(y),initial_relative_L2=norm(a-x)/norm(x),
            error_in_change_relative_L2=norm(change_error)/den,rms_error=norm(error)/np.sqrt(v.sum()),
            maximum_error=float(abs(error).max()),axis_rms_error=float(np.sqrt(np.sum(v[0]*error[0]**2)/v[0].sum())))
    return result
