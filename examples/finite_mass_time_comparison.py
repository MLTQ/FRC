"""Physical differences between archived trajectories on one fixed mesh."""
import numpy as np
from plasma_bench.finite_mass_state import evaluate
from .finite_mass_audit import load,state_from
from .finite_mass_evolution import case_state,protocol
from .current_thermal_norms import vectors,scalar


def difference(directory,first,second):
    cases={case['name']:case for case in protocol()}
    ref,mesh,initial,par,start,source=case_state(cases[first])
    a,aa=load(directory,first);b,bb=load(directory,second)
    assert a['completed'] and b['completed'] and a['actual_time_s']==b['actual_time_s']
    assert a['configuration']['nr']==b['configuration']['nr']
    a,b=state_from(aa,'final_'),state_from(bb,'final_')
    ea,eb=evaluate(mesh,a,par),evaluate(mesh,b,par)
    pa=ea.model.electron.pressure(a.electron_entropy);pb=eb.model.electron.pressure(b.electron_entropy)
    pi=(par.gamma_i-1)*(ea.components['ion']-eb.components['ion'])
    sums=dict(magnetic=0.,current=0.,electron_pressure=0.,ion_pressure=float(np.sum(mesh.volumes*pi*pi)))
    nodes,weights=np.polynomial.legendre.leggauss(10)
    for x,wx in zip(nodes,weights):
        r=mesh.r[:,None]+.5*mesh.dr*x
        for z,wz in zip(nodes,weights):
            az=.5*(1+z)
            # Linear field evaluation of differences avoids subtracting a bias.
            db,dj,ax=vectors(mesh,a.psi-b.psi,a.g-b.g,ea.current-eb.current,r,az)
            pe=scalar(pa-pb,ax,az);weight=.5*np.pi*wx*wz*r*mesh.dr*mesh.dz
            for key,value in (('magnetic',np.sum(db*db,axis=-1)),('current',np.sum(dj*dj,axis=-1)),('electron_pressure',pe*pe)):
                sums[key]+=float(np.sum(weight*value))
    return {key:float(np.sqrt(value/mesh.volumes.sum())) for key,value in sums.items()}


def measure(directory):
    result={}
    for family,names in (('manufactured',['moving_time_coarse','moving_16','moving_time_half','moving_time_quarter']),
                         ('fast_wave',['wave_time_coarse','wave_16','wave_time_half'])):
        pairs=[dict(first=a,second=b,rms_difference=difference(directory,a,b)) for a,b in zip(names[:-1],names[1:])]
        ratios=[{key:a['rms_difference'][key]/b['rms_difference'][key] if b['rms_difference'][key]>0 else None
                 for key in a['rms_difference']} for a,b in zip(pairs[:-1],pairs[1:])]
        result[family]=dict(pairs=pairs,successive_difference_ratios=ratios)
    return result
