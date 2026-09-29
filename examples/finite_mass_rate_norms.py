"""Independent parent conservative and species-pressure rate error norms."""
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.finite_mass_rhs import pressure_rates
from .finite_mass_manufactured import point_conserved,cell_average


def measure(mesh,ref,state,parameters,rate,time,*,family,order=8):
    if family=='manufactured':
        exact=cell_average(mesh,lambda r,z:point_conserved(ref,r,z,time,rate=True),order)
        piref=cell_average(mesh,lambda r,z:ref.derivatives(r,z,time)[2]['pi'],order)
        peref=ref.derivatives(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'),time)[2]['pe']
    else:
        def conservative(r,z):
            f=ref.derivatives(r,z,time)[2]
            energy=f['electron_pressure']/(ref.gamma-1)+f['ion_pressure']/(ref.gamma_i-1)+ref.bias*f['magnetic'][...,2]/MU0
            return np.stack((f['density'],ref.density*f['velocity'][...,0],ref.density*f['velocity'][...,2],energy,ref.density*r*f['velocity'][...,1]),axis=-1).real
        exact=cell_average(mesh,conservative,order)
        piref=cell_average(mesh,lambda r,z:ref.derivatives(r,z,time)[2]['ion_pressure'].real,order)
        peref=ref.derivatives(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'),time)[2]['electron_pressure'].real
    pedot,pidot=pressure_rates(state,rate,parameters)
    values={key:(rate.state.conserved[...,i],exact[...,i],mesh.volumes) for i,key in enumerate(('mass','radial_momentum','axial_momentum','total_energy','angular_momentum'))}
    values['electron_pressure']=(pedot,peref,rate.evaluation.model.mapping.volumes)
    values['ion_pressure']=(pidot,piref,mesh.volumes)
    metrics={}
    for key,(value,target,volume) in values.items():
        error=value-target
        metrics[key]=dict(relative_L2=float(np.sqrt(np.sum(volume*error**2)/np.sum(volume*target**2))),
            axis_relative_L2=float(np.sqrt(np.sum(volume[0]*error[0]**2)/np.sum(volume[0]*target[0]**2))),
            rms_error=float(np.sqrt(np.sum(volume*error**2)/np.sum(volume))),maximum_error=float(abs(error).max()))
    return metrics,dict(exact_conservative_rate=exact,exact_electron_pressure_rate=peref,exact_ion_pressure_rate=piref,
                        electron_pressure_rate=pedot,ion_pressure_rate=pidot)
