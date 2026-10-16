"""Evaluate candidate collision moments on unchanged completed pulse states."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E, electron_mass as ME, epsilon_0 as EPS0
from plasma_bench.material_experiment import build
from plasma_bench.maxwellian_pair_collisions import pair_moments
from plasma_bench.coulomb_transport import coulomb_coefficients
from .material_coil_plot import state_at
from .implicit_plasma_records import provenance, freeze, sha


def density_log_gradient(model, state, response):
    """Physical gradient of the declared n0/det(F) profile inside each Q1 cell."""
    f=model.fluid;mesh=model.space.mesh;corners=state.fluid.coordinates[f.ids]
    a=np.einsum('qnb,cna->cqab',f.reference_gradient,corners)/np.array([mesh.dr,mesh.dz])
    inverse=np.linalg.inv(a)
    mixed=(corners[:,0]-corners[:,1]+corners[:,2]-corners[:,3])/(mesh.dr*mesh.dz)
    r=model.space.reference_r.reshape(f.particles.shape);z=model.space.reference_z.reshape(f.particles.shape)
    radius=response['fields']['points'][...,0]
    grad_log_j=np.stack((a[...,0,0]/radius-1/r+np.einsum('cqa,ca->cq',inverse[...,1,:],mixed),
                         a[...,0,1]/radius+np.einsum('cqa,ca->cq',inverse[...,0,:],mixed)),axis=-1)
    grad_log_n0=np.stack((-3*r/(mesh.radius**2-r*r),-3*z/(mesh.half_length**2-z*z)),axis=-1)
    return np.linalg.solve(a.swapaxes(-1,-2),(grad_log_n0-grad_log_j)[...,None])[...,0].reshape(-1,2)


def regime_content(model, state, response, moments):
    """Content of conditional transport review regions, without changing state."""
    n=response['fields']['n'].ravel();te,ti=response['fields']['temperature'].reshape(-1,2).T
    volume=response['volume'];field=response['magnetic_field'];strength=np.linalg.norm(field,axis=1)
    gradient=density_log_gradient(model,state,response)
    inverse_length=np.linalg.norm(gradient,axis=1)
    parallel_inverse=np.divide(abs(np.sum(field[:,[0,2]]*gradient,axis=1)),strength,
                               out=np.zeros_like(strength),where=strength>0)
    coefficients=[]
    for density,electron,ion in zip(n,te,ti):
        try:coefficients.append(coulomb_coefficients(density,electron,ion,model.fluid.ion_mass))
        except ValueError:coefficients.append(None)
    available=np.array([value is not None for value in coefficients])
    column=lambda key:np.array([value[key] if value is not None else np.nan for value in coefficients])
    log_ei=column('coulomb_log_ei');log_ii=column('coulomb_log_ii')
    mfp_e=column('electron_thermal_mean_free_path_m');mfp_i=column('ion_thermal_mean_free_path_m')
    hall_e=E*strength/ME/column('electron_ion_momentum_rate_s')
    debye=np.sqrt(EPS0*E*te/(n*E**2))
    gyro=np.divide(np.sqrt(E*te/ME),E*strength/ME,out=np.full_like(strength,np.inf),where=strength>0)
    masks=dict(thermal_coefficients_unavailable=~available,
        small_thermal_coulomb_log=available&((log_ei<5)|(log_ii<5)),
        electron_unmagnetized_density_knudsen=available&(hall_e<1)&(mfp_e*inverse_length>=.1),
        electron_magnetized_parallel_density_knudsen=available&(hall_e>=1)&(mfp_e*parallel_inverse>=.1),
        ion_density_knudsen=available&(mfp_i*inverse_length>=.1),
        electron_drift_over_thermal=np.linalg.norm(response['drift'],axis=1)>=np.sqrt(E*te/ME),
        electron_gyroradius_near_debye=gyro<=10*debye)
    weights=dict(volume=volume,particles=volume*n,electron_thermal_energy=volume*n*te,
                 candidate_drag_power=volume*moments['drag_power'])
    content={name:dict(fractions={key:float(weight[mask].sum()/weight.sum()) if weight.sum() else None
                                 for key,weight in weights.items()},
                       candidate_drag_power_w=float(weights['candidate_drag_power'][mask].sum()))
             for name,mask in masks.items()}
    ranges={name:([float(values[available].min()),float(values[available].max())] if available.any() else None)
            for name,values in [('thermal_log_ei',log_ei),('thermal_log_ii',log_ii),
                                ('electron_mfp_over_density_length',mfp_e*inverse_length),
                                ('electron_mfp_over_parallel_density_length',mfp_e*parallel_inverse),
                                ('electron_hall_parameter',hall_e)]}
    return dict(content=content,ranges=ranges,
        scope='Conditional thermal NRL coefficients at actual positive plasma samples; zero-drift rates are not a finite-drift cutoff model. Density gradients are within-cell analytic derivatives; entropy/field gradients, interface jumps and nonlocal distributions are not resolved by these indicators. Thresholds are review criteria, never a validity certificate; regions overlap.')


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
            max_power_pair_residual_w_m3=float(abs(moments['heat_to_a']+moments['heat_to_b']-moments['drag_power']).max()),
            regime=regime_content(model,state,a,moments)))
    t=np.array([r['time_s'] for r in rows]);result=dict(source=str(source),inputs=inputs,
        diagnostic_source_sha256=identity['source_sha256'],pulse_source_sha256=data['source_sha256'],
        declared_coulomb_log=10.,prescribed_resistivity_ohm_m=eta,records=rows,
        scope='Postprocessing only: actual fixed states from the prescribed-eta pulse, with hypothetical isotropic Maxwellian lnLambda=10 collision moments. No new trajectory, physical cutoff selection, kinetic validity or predicted total heat under changed transport.')
    result['path_integrals_j']={key:float(np.trapezoid([r[key] for r in rows],t)) for key in
        ('prescribed_joule_power_w','candidate_drag_power_w','linear_maxwellian_drag_power_w','candidate_electron_heat_power_w','candidate_ion_heat_power_w')}
    result['review_region_path_energy_fraction']={name:float(np.trapezoid(
        [r['regime']['content'][name]['candidate_drag_power_w'] for r in rows],t)/result['path_integrals_j']['candidate_drag_power_w'])
        if result['path_integrals_j']['candidate_drag_power_w'] else None for name in rows[0]['regime']['content']}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
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
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),constrained_layout=True)
    categories=[('small_thermal_coulomb_log','Thermal Coulomb log <5'),
        ('electron_unmagnetized_density_knudsen','Electron density Kn ≥0.1, weak B'),
        ('electron_magnetized_parallel_density_knudsen','Electron parallel density Kn ≥0.1, strong B'),
        ('ion_density_knudsen','Ion density Kn ≥0.1'),('electron_drift_over_thermal','Electron drift ≥thermal speed'),
        ('thermal_coefficients_unavailable','Thermal coefficients unavailable')]
    for axis,weight,label in zip(axes,('particles','electron_thermal_energy','candidate_drag_power'),
                                 ('Particle inventory','Electron thermal energy','Candidate drag power')):
        for key,name in categories:
            values=[r['regime']['content'][key]['fractions'][weight] for r in rows]
            axis.plot(t*1e6,[None if v is None else 100*v for v in values],label=name)
        axis.set(xlabel='Time [µs]',ylabel='Fraction [%]',title=label,ylim=(0,105));axis.grid(alpha=.2)
    axes[-1].legend(fontsize=7,loc='center left',bbox_to_anchor=(1,0.5))
    fig.suptitle('Conditional transport review across the actual plasma\nOverlapping regions; density scales only; thermal estimates do not establish finite-drift validity')
    fig.savefig(out/'collision_regime_content.png',dpi=170);fig.savefig(out/'collision_regime_content.pdf');plt.close(fig)
    if provenance(['examples.material_collision_diagnostic'])!=identity:raise RuntimeError('Diagnostic source changed')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('Collision sensitivity complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source');parser.add_argument('destination')
    args=parser.parse_args();diagnose(args.source,args.destination)
