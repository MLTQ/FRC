"""Full existing atomic histories with one explicit private cross-section change."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV, electron_mass as ME
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.cross_section_ionization import ArgonGroundCrossSection
from plasma_bench.lxcat_cross_sections import CrossSection, MaxwellianIntegral
from .atomic_network import argon_models, evolve
from .implicit_plasma_records import provenance, freeze


def run(destination, neutral_path, ion_path, ionization_path, bsr_path):
    out = Path(destination)
    if out.exists(): raise ValueError('Fresh destination required')
    controls, inputs = argon_models(neutral_path, ion_path, ionization_path)
    control = controls['joint']
    replacement = ArgonGroundCrossSection(control.ionization, bsr_path)
    mappings = [(int(control.charge[origin]), int(np.flatnonzero(control.blocks[control.charge[origin]] == origin)[0]), row)
                for origin, row in control.resolved_ionization.items()]
    model = AtomicNetwork(replacement, control.recombination, control.levels, resolved_ionization=mappings)
    identity = provenance(['examples.bsr_atomic_source']); out.mkdir(parents=True); freeze(out, identity)
    started = time.monotonic(); histories = {}
    report = dict(terminal=False, source_sha256=identity['source_sha256'], inputs=inputs,
        bsr_input={k:replacement.data[k] for k in ('input_path','input_sha256')}, cases=[], comparisons=[],
        scope='Homogeneous level/charge/electron-energy source. Only ground-neutral Ar ionization changes to a finite-domain BSR Maxwellian integral. No moving or current coupling, excited-parent completion, photon escape or formation claim.',
        input_policy='Private raw tables and reconstructed coefficient arrays remain outside repository. Only code, input identities and actual simulation histories are archived.')
    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    # Independent integral for sigma(E)=a+bE on a finite interval. This evaluates
    # the physical Maxwellian normalization without using the table quadrature.
    energies = np.array([2.,3.,5.,11.,25.,80.]); a=2e-20; b=1e-22
    synthetic = CrossSection('IONIZATION','manufactured',np.array([2.]),(),energies,a+b*energies)
    quadrature = []
    for order in (4,8,16):
        integrated = MaxwellianIntegral(synthetic,order); errors=[]
        for temperature in (1.,2.,5.,10.,20.,40.):
            lo,hi=energies[[0,-1]]/temperature
            f1=lambda x:(x+1)*np.exp(-x)
            f2=lambda x:(x*x+2*x+2)*np.exp(-x)
            exact=np.sqrt(8*EV*temperature/(np.pi*ME))*(a*(f1(lo)-f1(hi))+b*temperature*(f2(lo)-f2(hi)))
            errors.append(abs(integrated(temperature)/exact-1))
        quadrature.append(dict(order=order,maximum_relative_error=float(max(errors))))
    section = replacement.rate.section; lower_order=MaxwellianIntegral(section,8)
    report['quadrature_reference'] = quadrature
    report['physical_integral_order8_16_max_relative_difference'] = max(
        abs(lower_order(t)/replacement.rate(t)-1) for t in (1.,2.,5.,10.,20.,40.))
    report['table_domain_ev'] = section.energy_ev[[0,-1]].tolist()
    report['source_threshold_ev'] = float(section.parameters[0])
    report['physical_ionization_energy_ev'] = float(replacement.ionization_energy[0]/EV)
    report['ground_rate_ratio_bsr_to_control_at_10ev'] = replacement.rate(10.)/control.ionization.coefficients(10*EV)[0][0]
    report['range_note'] = 'No validated Maxwellian temperature interval is supplied: replacement events are flagged outside/unspecified fit range. Above-table electron flux fraction is not a missing-rate error bound.'
    save()
    density = 8.961429541213399e19
    for case,seed,deposit in (('singly_ionized',1.,0.),('partial_seed',.1,0.),('heated_seed',.1,20.)):
        fraction=np.zeros(19);fraction[0]=1-seed;fraction[1]=seed
        for name,source,step_counts in (('control',control,(1400,)),('bsr',model,(350,700,1400))):
            for steps in step_counts:
                summary,history=evolve(source,source.initial(fraction,10.),density,3.5e-6,steps,deposit)
                key=f'{case}_{name}_{steps}';histories[key]=history
                np.savez_compressed(out/f'{key}.npz',**history)
                max_temperature=float(history['temperature_ev'].max())
                summary.update(case=case,model=name,steps=steps,nucleus_density_m3=density,
                    initial_singly_ionized_fraction=seed,initial_temperature_ev=10.,deposit_ev_per_nucleus=deposit,
                    maximum_temperature_ev=max_temperature,
                    maximum_electron_flux_fraction_above_bsr_table=replacement.rate.electron_flux_fraction_above_table(max_temperature))
                report['cases'].append(summary);save()
                print(case,name,steps,summary['success'],summary['final_mean_charge'],summary['final_temperature_ev'],flush=True)
        def compare(first,second,stride,label):
            x=histories[first];y=histories[second]
            if len(x['time_s']) != len(y['time_s'][::stride]) or not np.allclose(x['time_s'],y['time_s'][::stride],rtol=0,atol=1e-20):
                return
            report['comparisons'].append(dict(case=case,comparison=label,
                maximum_population_l1_difference=float(np.abs(x['states'][:,:model.count]-y['states'][::stride,:model.count]).sum(axis=1).max()),
                maximum_temperature_difference_ev=float(abs(x['temperature_ev']-y['temperature_ev'][::stride]).max()),
                maximum_charge_difference=float(abs(x['mean_charge']-y['mean_charge'][::stride]).max()),
                final_temperature_difference_ev=float(y['temperature_ev'][-1]-x['temperature_ev'][-1]),
                final_mean_charge_difference=float(y['mean_charge'][-1]-x['mean_charge'][-1])))
        completed={f"{c['case']}_{c['model']}_{c['steps']}" for c in report['cases'] if c['success']}
        for first,second,stride,label in ((f'{case}_bsr_350',f'{case}_bsr_700',2,'bsr_10_to_5ns'),
                (f'{case}_bsr_700',f'{case}_bsr_1400',2,'bsr_5_to_2.5ns'),
                (f'{case}_control_1400',f'{case}_bsr_1400',1,'ground_rate_source_only')):
            if first in completed and second in completed: compare(first,second,stride,label)
        save()
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for column,case in enumerate(('singly_ionized','partial_seed','heated_seed')):
        for name,label in (('control','Existing ground-Ar fit'),('bsr','BSR ground-Ar cross section')):
            data=histories[f'{case}_{name}_1400'];t=data['time_s']*1e6
            axes[0,column].plot(t,data['temperature_ev'],label=label)
            axes[1,column].plot(t,data['mean_charge'],label=label)
        axes[0,column].set_title(case.replace('_',' '))
        for row,label in enumerate(('Electron temperature [eV]','Electrons per nucleus')):
            axes[row,column].set(xlabel='Time [µs]',ylabel=label);axes[row,column].grid(alpha=.2)
    axes[0,0].legend(fontsize=9)
    fig.suptitle('One changed atomic input: ground-state argon ionization\nHomogeneous source comparison; produced photons are not established escape')
    fig.savefig(out/'bsr_atomic_source.png',dpi=180);fig.savefig(out/'bsr_atomic_source.pdf');plt.close(fig)
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save()
    return report['success']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');p.add_argument('neutral_path')
    p.add_argument('ion_path');p.add_argument('ionization_path');p.add_argument('bsr_path')
    args=p.parse_args()
    raise SystemExit(0 if run(args.destination,args.neutral_path,args.ion_path,args.ionization_path,args.bsr_path) else 1)
