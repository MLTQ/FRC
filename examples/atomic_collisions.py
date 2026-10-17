"""Complete reactive collision histories and archived independent collision limit."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, electron_mass as ME
from plasma_bench.atomic_collisions import CollisionalAtomicMoments
from plasma_bench.elastic_atomic import ElasticAtomicMoments
from plasma_bench.lxcat_cross_sections import read_cross_sections
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization, evolve
from .atomic_network import argon_models
from .implicit_plasma_records import provenance, freeze


def plot_histories(out, plotted):
    """Show all current components needed to distinguish damping from redirection."""
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for name,h in plotted:
        label=name+(' (incomplete)' if len(h['states'])!=1401 else '')
        values=(h['electrons'],h['electron_temperature_ev'],np.linalg.norm(h['current_a_m2'],axis=1),
                h['current_a_m2'][:,1],h['current_a_m2'][:,2],h['energy_budget_ev']-h['energy_budget_ev'][0])
        for ax,y in zip(axes.flat,values):
            ax.plot(h['time_s']*1e6,y,label=label)
    for ax,label in zip(axes.flat,('Electrons / nucleus','Electron temperature [eV]','Current magnitude [A/m²]',
            'Azimuthal current [A/m²]','Axial current [A/m²]','Energy residual [eV/nucleus]')):
        ax.set(xlabel='Time [µs]',ylabel=label);ax.grid(alpha=.2)
    axes[0,2].set_yscale('log')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Heated argon: reacting populations with paired collision force and heat\nHomogeneous source; prescribed lnΛ=10; neutral and ion flows initially differ',fontsize=12)
    fig.savefig(out/'reactive_collisions.png',dpi=170);fig.savefig(out/'reactive_collisions.pdf');plt.close(fig)


def run(destination, collision_reference, reaction_reference, neutral_path, ion_path, ionization_path,
        neutral_elastic_file=None, charged_reference=None):
    out=Path(destination); collision_reference=Path(collision_reference); reaction_reference=Path(reaction_reference)
    if out.exists():raise ValueError('Fresh output directory required')
    old_collision=json.loads((collision_reference/'summary.json').read_text())
    old_reaction=json.loads((reaction_reference/'summary.json').read_text())
    if not old_collision['terminal'] or not old_reaction['terminal']:
        raise ValueError('Terminal archived reference runs required')
    models,inputs=argon_models(neutral_path,ion_path,ionization_path)
    elastic_section = None
    if neutral_elastic_file is not None:
        if charged_reference is None:
            raise ValueError('Neutral-collision comparison requires the complete charged-only reference')
        charged_reference = Path(charged_reference)
        charged_report = json.loads((charged_reference/'summary.json').read_text())
        if not charged_report.get('terminal') or not charged_report.get('success'):
            raise ValueError('Complete successful charged-only reference required')
        data = read_cross_sections(neutral_elastic_file)
        sections = [s for s in data['sections'] if s.kind == 'ELASTIC' and s.target == 'Ar']
        if len(sections) != 1:
            raise ValueError('Exactly one explicit Ar ELASTIC momentum-transfer table required')
        elastic_section = sections[0]
        inputs['neutral_elastic'] = dict(path=data['input_path'], sha256=data['input_sha256'],
            provider='LXCat BSR-500, O. Zatsarinny and K. Bartschat; https://www.lxcat.net/BSR',
            energy_bounds_ev=elastic_section.energy_ev[[0, -1]].tolist())
    identity=provenance(['examples.atomic_collisions']);out.mkdir(parents=True);freeze(out,identity)
    report=dict(terminal=False,source_sha256=identity['source_sha256'],inputs=inputs,cases=[],
        collision_reference=str(collision_reference.resolve()),collision_reference_source=old_collision['source_sha256'],
        reaction_reference=str(reaction_reference.resolve()),reaction_reference_source=old_reaction['source_sha256'],
        neutral_elastic_enabled=elastic_section is not None,
        scope='Homogeneous joint atomic reactions plus all charged Maxwellian collision pairs. Explicit lnLambda=10. Ground-neutral/electron elastic scattering only when explicitly enabled; other neutral channels, charge exchange, spatial flux, electromagnetic evolution and photon escape remain absent.')
    if elastic_section is not None:
        report.update(charged_reference=str(charged_reference.resolve()),
            charged_reference_source=charged_report['source_sha256'],
            collision_limits_reused='Original complete collision-only limits are archived in charged_reference; no rerun of unaffected limits.')
    started=time.monotonic();histories={}
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def record(name,model,initial,density,duration,steps,deposit=0.,reference=None):
        if isinstance(model, ElasticAtomicMoments):
            model.maximum_elastic_outside_probability = 0.
        summary,history=evolve(model,initial,density,duration,steps,deposit)
        summary.update(name=name,coulomb_log=model.coulomb_log,
            reference_velocity_m_s=model.reference_velocity.tolist(),
            maximum_cold_charged_populations=int(history['cold_charged_populations'].max()),
            final_mean_ion_temperature_ev=float(history['mean_ion_temperature_ev'][-1]),
            final_mean_neutral_temperature_ev=float(history['mean_neutral_temperature_ev'][-1]))
        if isinstance(model, ElasticAtomicMoments):
            neutral = np.flatnonzero(model.network.charge == 0)
            missing = neutral[neutral != model.neutral_index]
            total = history['states'][:, neutral].sum(axis=1)
            excluded = history['states'][:, missing].sum(axis=1)
            history['excited_neutral_without_elastic_per_nucleus'] = excluded
            summary.update(maximum_elastic_probability_outside_table=model.maximum_elastic_outside_probability,
                maximum_excited_neutral_without_elastic_per_nucleus=float(excluded.max()),
                maximum_neutral_fraction_without_elastic=float(np.divide(excluded,total,out=np.zeros_like(total),where=total>0).max()))
        if reference is not None:
            primitive=[]
            for state in history['states']:
                f=model.fields(state);v=f['velocity']+model.reference_velocity
                primitive.append(np.r_[v[-1],v[1],f['temperature_ev'][-1],f['temperature_ev'][1]])
            primitive=np.asarray(primitive);expected=reference[:len(primitive)]
            history.update(primitive=primitive,independent_reference=expected)
            delta=(primitive[:,:3]-primitive[:,3:6])-(expected[:,:3]-expected[:,3:6])
            summary['independent_comparison']=dict(
                maximum_relative_drift_error=float(np.linalg.norm(delta,axis=1).max()/np.linalg.norm(expected[0,:3]-expected[0,3:6])),
                maximum_temperature_errors_ev=abs(primitive[:,6:]-expected[:,6:]).max(axis=0).tolist())
        np.savez_compressed(out/f'{name}.npz',**history)
        histories[name]=history;report['cases'].append(summary);save()
        print('Atomic collisions',name,summary['success'],summary.get('failure'),summary.get('independent_comparison'),flush=True)
    save()
    # Match the existing physical electron/argon collision problem exactly.
    # Its ion mass is 39.948 u; add me to the neutral mass so m_ion is unchanged.
    rr=RadiativeRecombination(dict(nuclear_charge=1,transitions=[dict(initial_charge=1,kind='power_law',
        A_cm3_s=1e-12,reference_temperature_ev=2.,exponent=0.,minimum_temperature_ev=None,maximum_temperature_ev=None)]))
    network=AtomicNetwork(ConstantIonization(),rr,ionization_enabled=False,recombination_enabled=False,radiation=False)
    case=old_collision['cases']['electron_argon'];ma,mb=case['masses_kg']
    if ma!=ME or case['density_m3'][0]!=case['density_m3'][1]:raise ValueError('Reference requires physical electrons and equal singly charged densities')
    model=CollisionalAtomicMoments(network,mb+ME,[200.,-100.,300.],coulomb_log=case['coulomb_log'])
    for steps in (() if elastic_section is not None else (100,200,400)):
        with np.load(collision_reference/f'electron_argon_{steps}.npz') as data:reference=data['reference'].copy()
        initial=model.initial(network.initial([0.,1.],reference[0,6]),reference[0,7],
            np.vstack((np.zeros(3),reference[0,3:6],reference[0,:3])))
        record(f'collision_limit_{steps}',model,initial,case['density_m3'][0],case['duration_s'],steps,reference=reference)
    network=models['joint']
    model=(CollisionalAtomicMoments(network,39.948*atomic_mass,[0.,0.,0.],coulomb_log=10.)
           if elastic_section is None else
           ElasticAtomicMoments(network,39.948*atomic_mass,[0.,0.,0.],coulomb_log=10.,elastic_section=elastic_section))
    velocity=np.zeros((network.count+1,3));velocity[:-1,2]=np.where(network.charge==0,1000.,-1000.)
    velocity[-1]=[0.,10000.,0.]
    initial=model.initial(network.initial(np.r_[.9,.1,np.zeros(17)],10.),.1,velocity)
    for steps in (700,1400):record(f'heated_seed_{steps}',model,initial,8.961429541213399e19,3.5e-6,steps,20.)
    coarse,fine=histories['heated_seed_700'],histories['heated_seed_1400']
    if len(coarse['states'])==701 and len(fine['states'])==1401:
        report['time_halving']=dict(maximum_population_l1=float(abs(coarse['states'][:,:network.count]-fine['states'][::2,:network.count]).sum(axis=1).max()),
            maximum_electron_temperature_difference_ev=float(abs(coarse['electron_temperature_ev']-fine['electron_temperature_ev'][::2]).max()),
            maximum_current_difference_a_m2=float(np.linalg.norm(coarse['current_a_m2']-fine['current_a_m2'][::2],axis=1).max()))
    with np.load(reaction_reference/'heated_seed_1400.npz') as data:previous={k:data[k].copy() for k in data.files}
    report['reaction_only_final']=dict(electrons=float(previous['electrons'][-1]),
        temperature_ev=float(previous['electron_temperature_ev'][-1]),current_a_m2=previous['current_a_m2'][-1].tolist())
    plotted = [('Reactions only',previous)]
    if elastic_section is not None:
        with np.load(charged_reference/'heated_seed_1400.npz') as data:
            charged = {k:data[k].copy() for k in data.files}
        plotted.append(('Reactions + charged collisions', charged))
        if len(fine['states']) == 1401 and np.array_equal(fine['time_s'], charged['time_s']):
            report['charged_only_comparison'] = dict(
                maximum_population_l1=float(abs(fine['states'][:,:network.count]-charged['states'][:,:network.count]).sum(axis=1).max()),
                maximum_electron_temperature_difference_ev=float(abs(fine['electron_temperature_ev']-charged['electron_temperature_ev']).max()),
                maximum_current_difference_a_m2=float(np.linalg.norm(fine['current_a_m2']-charged['current_a_m2'],axis=1).max()),
                final_current_difference_a_m2=(fine['current_a_m2'][-1]-charged['current_a_m2'][-1]).tolist(),
                final_neutral_temperature_difference_ev=float(fine['mean_neutral_temperature_ev'][-1]-charged['mean_neutral_temperature_ev'][-1]))
    plotted.append(('Reactions + charged + ground-neutral collisions' if elastic_section is not None else 'Reactions + collisions', fine))
    plot_histories(out,plotted)
    if provenance(['examples.atomic_collisions'])!=identity:raise RuntimeError('Source changed during integration')
    report['terminal']=True;report['success']=all(c['success'] for c in report['cases']);save()
    print('Atomic collisions terminal',report['success'],out,flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('destination','collision_reference','reaction_reference','neutral_input','ion_input','ionization_input'):parser.add_argument(name)
    parser.add_argument('--neutral-elastic-file');parser.add_argument('--charged-reference')
    a=parser.parse_args();report=run(a.destination,a.collision_reference,a.reaction_reference,a.neutral_input,a.ion_input,a.ionization_input,
                                   a.neutral_elastic_file,a.charged_reference)
    raise SystemExit(0 if report['success'] else 1)
