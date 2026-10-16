"""Pressure-driven neutral slab expansion against the exact vacuum rarefaction."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV
from plasma_bench.atomic_transport import AtomicTransport
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out=Path(destination)
    if out.exists():raise ValueError('Fresh expansion output directory required')
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    network=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=False,recombination_enabled=False,radiation=False)
    model=FixedFrameAtomicMoments(network,39.948*atomic_mass,[0.,0.,0.])
    gamma=5/3;length=.2;half_width=.025;temperature=.1
    sound=np.sqrt(gamma*EV*temperature/model.mass[0]);duration=half_width/(4*sound)
    identity=provenance(['examples.atomic_expansion']);out.mkdir(parents=True);freeze(out,identity)
    report=dict(terminal=False,source_sha256=identity['source_sha256'],cases=[],
        length_m=length,initial_half_width_m=half_width,initial_temperature_ev=temperature,duration_s=duration,
        scope='Neutral Euler vacuum-rarefaction component; no electrons, fields, collisions or formation. Numerical density leakage beyond the exact front is measured, never called an exact material boundary.')
    started=time.monotonic();finals=[]
    def save():
        report['elapsed_seconds']=time.monotonic()-started;(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def reference(transport,t):
        result=np.zeros((transport.cells,model.size));nodes,weights=np.polynomial.legendre.leggauss(8)
        if t==0:breaks=[-half_width,half_width]
        else:breaks=[-half_width-2*sound*t/(gamma-1),-half_width+sound*t,half_width-sound*t,half_width+2*sound*t/(gamma-1)]
        for cell,centre in enumerate(transport.x-length/2):
            lo,hi=centre-transport.dx/2,centre+transport.dx/2
            cuts=np.r_[lo,[x for x in breaks if lo<x<hi],hi]
            for left,right in zip(cuts[:-1],cuts[1:]):
                x=.5*(left+right)+.5*(right-left)*nodes
                if t==0:
                    density=(abs(x)<half_width).astype(float);pressure=temperature*density;velocity=np.zeros_like(x)
                else:
                    xi=(abs(x)-half_width)/t
                    core=xi<=-sound;fan=(xi>-sound)&(xi<2*sound/(gamma-1))
                    ratio=np.zeros_like(x);ratio[core]=1.;ratio[fan]=2/(gamma+1)*(1-(gamma-1)*xi[fan]/(2*sound))
                    density=ratio**(2/(gamma-1));pressure=temperature*ratio**(2*gamma/(gamma-1))
                    velocity=np.zeros_like(x);velocity[fan]=np.sign(x[fan])*2/(gamma+1)*(sound+xi[fan])
                w=.5*(right-left)*weights/transport.dx
                result[cell,0]+=w@density
                result[cell,model.momentum[0,2]]+=w@(density*velocity/model.speed_scale)
                result[cell,model.heavy_energy[0]]+=w@(pressure/(gamma-1)+.5*(model.mass[0]/EV)*density*velocity**2)
        return result
    save()
    for cells in (64,128,256):
        transport=AtomicTransport(model,cells,length,1e19);state=reference(transport,0.)
        t=0.;times=[0.];states=[state.copy()];failure=None
        while t<duration:
            dt=min(transport.timestep(state),duration-t)
            try:candidate=transport.transport_step(state,dt)
            except (ValueError,RuntimeError) as error:
                failure=dict(step=len(times),time_s=t,dt_s=dt,message=str(error));break
            state=candidate;t=min(t+dt,duration);times.append(t);states.append(state.copy())
        states=np.asarray(states);expected=reference(transport,t)
        totals=[transport.totals(s) for s in states];budgets={k:np.asarray([r[k] for r in totals]) for k in totals[0]}
        front=half_width+2*sound*t/(gamma-1);outside=abs(transport.x-length/2)-transport.dx/2>=front
        case=dict(cells=cells,steps=len(times)-1,success=failure is None,failure=failure,completed_time_s=t,
            final_density_mean_l1=float(abs(state[:,0]-expected[:,0]).mean()),
            final_momentum_mean_l1_scaled=float(abs(state[:,model.momentum[0,2]]-expected[:,model.momentum[0,2]]).mean()),
            final_translational_energy_mean_l1_ev=float(abs(state[:,model.heavy_energy[0]]-expected[:,model.heavy_energy[0]]).mean()),
            maximum_number_integral_error=float(abs(budgets['nuclei']-budgets['nuclei'][0]).max()),
            maximum_momentum_integral_error=float(np.linalg.norm(budgets['momentum_scaled']-budgets['momentum_scaled'][0],axis=1).max()),
            maximum_energy_integral_error_ev_m=float(abs(budgets['energy_ev']-budgets['energy_ev'][0]).max()),
            number_fraction_in_wholly_exterior_cells=float(transport.dx*state[outside,0].sum()/budgets['nuclei'][0]),
            maximum_density_in_wholly_exterior_cells=float(state[outside,0].max(initial=0.)),
            exactly_empty_final_cells=int(np.count_nonzero(state[:,0]==0)),
            maximum_electron_number=float(states[:,:,network.electrons].max()))
        np.savez_compressed(out/f'neutral_expansion_{cells}.npz',time_s=times,x_m=transport.x-length/2,states=states,analytic_final=expected,**budgets)
        report['cases'].append(case);finals.append((transport.x-length/2,state,expected));save();print('Neutral expansion',case,flush=True)
    fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    for cells,(x,state,expected) in zip((64,128,256),finals):
        for ax,index in zip(axes,(0,model.momentum[0,2],model.heavy_energy[0])):ax.plot(x,state[:,index],label=f'{cells} cells')
    for ax,index,label in zip(axes,(0,model.momentum[0,2],model.heavy_energy[0]),('Number / reference density','Scaled axial momentum','Translational energy [eV/reference nucleus]')):
        ax.plot(x,expected[:,index],'k--',label='Exact rarefaction cell average');ax.set(xlabel='Axial coordinate [m]',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Pressure-driven neutral expansion: Euler verification, with numerical front diffusion')
    fig.savefig(out/'neutral_expansion.png',dpi=170);fig.savefig(out/'neutral_expansion.pdf');plt.close(fig)
    if provenance(['examples.atomic_expansion'])!=identity:raise RuntimeError('Source changed during expansion')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save();return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');a=parser.parse_args()
    report=run(a.destination);raise SystemExit(0 if report['success'] else 1)
