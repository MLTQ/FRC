"""Independent continuous contact pressure/temperature and zero-current norms."""
import numpy as np
from scipy.constants import elementary_charge as E
from .current_thermal_norms import vectors,scalar


def measure(mesh,ref,state,par,evaluation,time,*,order=10):
    """Integrate the declared Q1 electron and cellwise ion representations."""
    model=evaluation.model
    pe=model.electron.pressure(state.electron_entropy)
    pi=(par.gamma_i-1)*evaluation.components['ion']
    ni=state.conserved[...,0]/par.masses.total
    keys=('electron_pressure','ion_pressure','electron_temperature_ev','ion_temperature_ev')
    error={k:0. for k in keys};norm=error.copy();maximum=error.copy()
    drift={k:0. for k in ('magnetic_t','current_a_m2')};peak=drift.copy()
    nodes,weights=np.polynomial.legendre.leggauss(order)
    current=model.current(model.geometry.pack(state.psi,state.g))
    for ar,wr in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*ar
        for az,wz in zip(nodes,weights):
            z=mesh.z[None]+.5*mesh.dz*az
            magnetic,j,x=vectors(mesh,state.psi,state.g,current,r,.5*(1+az))
            magnetic[...,2]-=ref.bias
            ne=scalar(model.electron.density,x,.5*(1+az))
            ep=scalar(pe,x,.5*(1+az))
            exact=ref.fields(z,time)
            actual=dict(electron_pressure=ep,ion_pressure=pi,
                electron_temperature_ev=ep/(E*ne),ion_temperature_ev=pi/(E*ni))
            target={k:exact[k] for k in ('electron_pressure','ion_pressure')}
            target.update(electron_temperature_ev=exact['electron_pressure']/(E*exact['number_density']),
                          ion_temperature_ev=exact['ion_pressure']/(E*exact['number_density']))
            weight=.5*np.pi*wr*wz*r*mesh.dr*mesh.dz
            for key in keys:
                defect=actual[key]-target[key]
                error[key]+=float(np.sum(weight*defect**2))
                norm[key]+=float(np.sum(weight*target[key]**2))
                maximum[key]=max(maximum[key],float(abs(defect).max()))
            for key,value in (('magnetic_t',magnetic),('current_a_m2',j)):
                squared=np.sum(value*value,axis=-1)
                drift[key]+=float(np.sum(weight*squared))
                peak[key]=max(peak[key],float(np.sqrt(squared.max())))
    volume=float(np.sum(mesh.volumes))
    return dict(fields={k:dict(relative_L2=float(np.sqrt(error[k]/norm[k])),
                rms_error=float(np.sqrt(error[k]/volume)),sampled_maximum_error=maximum[k]) for k in keys},
        zero_target_drift={k:dict(rms=float(np.sqrt(drift[k]/volume)),sampled_maximum=peak[k]) for k in drift})
