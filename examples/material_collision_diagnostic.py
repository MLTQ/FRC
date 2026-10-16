"""Evaluate candidate collision moments on unchanged completed pulse states."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E, electron_mass as ME
from plasma_bench.material_experiment import build
from plasma_bench.maxwellian_pair_collisions import pair_moments
from .material_coil_plot import state_at
from .implicit_plasma_records import provenance, freeze, sha


def diagnose(source, destination):
    source=Path(source);out=Path(destination);data=json.loads((source/'summary.json').read_text())
    if out.exists() or not data['terminal'] or not data['success']:
        raise ValueError('Complete pulse and fresh diagnostic output required')
    out.mkdir(parents=True);identity=provenance(['examples.material_collision_diagnostic']);freeze(out,identity)
    fluid,model,metadata=build(data['config']);eta=data['config']['plasma']['prescribed_resistivity_ohm_m']
    rows=[];inputs={str(source/'summary.json'):sha(source/'summary.json')}
    for recorded in [data['initial']]+data['records']:
        step=recorded.get('step',0);state=state_at(source,step);time=step*data['dt']
        raw=source/('initial.npz' if step==0 else f'step_{step:05}.npz');inputs[str(raw)]=sha(raw)
        a=model.response(state,rate=False,time=time);n=a['fields']['n'].ravel()
        te,ti=a['fields']['temperature'].reshape(-1,2).T;volume=a['volume'];w=a['drift']
        moments=pair_moments(n,n,te,ti,-w,ME,fluid.ion_mass,coulomb_log=10.)
        linear=pair_moments(n,n,te,ti,np.zeros_like(w),ME,fluid.ion_mass,coulomb_log=10.)
        drift2=np.sum(w*w,axis=1)
        prescribed=eta*(E*n)**2*drift2
        zero_drift=linear['drag_coefficient']*drift2
        frequency=moments['drag_coefficient']/n*(1/ME+1/fluid.ion_mass)
        rows.append(dict(step=step,time_s=time,prescribed_joule_power_w=float(volume@prescribed),
            candidate_drag_power_w=float(volume@moments['drag_power']),
            linear_maxwellian_drag_power_w=float(volume@zero_drift),
            candidate_electron_heat_power_w=float(volume@moments['heat_to_a']),
            candidate_ion_heat_power_w=float(volume@moments['heat_to_b']),
            maximum_relative_drift_relaxation_rate_s=float(frequency.max()),
            max_rate_times_original_dt=float(frequency.max()*data['dt']),
            finite_over_linear_power=float((volume@moments['drag_power'])/(volume@zero_drift)) if volume@zero_drift else None,
            max_power_pair_residual_w_m3=float(abs(moments['heat_to_a']+moments['heat_to_b']-moments['drag_power']).max())))
    t=np.array([r['time_s'] for r in rows]);result=dict(source=str(source),inputs=inputs,
        diagnostic_source_sha256=identity['source_sha256'],pulse_source_sha256=data['source_sha256'],
        declared_coulomb_log=10.,prescribed_resistivity_ohm_m=eta,records=rows,
        scope='Postprocessing only: actual fixed states from the prescribed-eta pulse, with hypothetical isotropic Maxwellian lnLambda=10 collision moments. No new trajectory, physical cutoff selection, kinetic validity or predicted total heat under changed transport.')
    result['path_integrals_j']={key:float(np.trapezoid([r[key] for r in rows],t)) for key in
        ('prescribed_joule_power_w','candidate_drag_power_w','linear_maxwellian_drag_power_w','candidate_electron_heat_power_w','candidate_ion_heat_power_w')}
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    for key,label in [('prescribed_joule_power_w','Original prescribed eta'),('linear_maxwellian_drag_power_w','Small-drift Maxwellian'),('candidate_drag_power_w','Finite-drift Maxwellian')]:
        axes[0].plot(t*1e6,np.array([r[key] for r in rows])/1e3,label=label)
    axes[0].set(ylabel='Instantaneous drag power [kW]',title='Same saved states; different collision laws')
    for key,label in [('candidate_electron_heat_power_w','Electron heat'),('candidate_ion_heat_power_w','Ion heat')]:
        axes[1].plot(t*1e6,np.array([r[key] for r in rows])/1e3,label=label)
    axes[1].set(ylabel='Hypothetical thermal power [kW]',title='Paired drag and temperature exchange')
    axes[2].plot(t*1e6,[r['max_rate_times_original_dt'] for r in rows],label='Maximum local rate × dt')
    axes[2].set(ylabel='Relative drift relaxation rate × original dt',title='Local collisional timescale screening')
    for ax in axes:ax.set_xlabel('Time [µs]');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Collision sensitivity on the completed long-array trajectory\nDeclared lnΛ=10 Maxwellians; no transport-modified evolution or validity claim')
    fig.savefig(out/'collision_sensitivity.png',dpi=170);fig.savefig(out/'collision_sensitivity.pdf');plt.close(fig)
    if provenance(['examples.material_collision_diagnostic'])!=identity:raise RuntimeError('Diagnostic source changed')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('Collision sensitivity complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source');parser.add_argument('destination')
    args=parser.parse_args();diagnose(args.source,args.destination)
