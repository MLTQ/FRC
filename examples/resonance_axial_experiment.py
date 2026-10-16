"""Argon detector decay with finite ends, integrated axial source and collisions."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as KB, c
from scipy.special import roots_legendre
from plasma_bench.line_radiation import RadiativeLevels
from plasma_bench.resonance_axial import axial_transport, resonant_decorrelation_rate
from .resonance_experiment import summarize_batches, vacuum_reference
from .implicit_plasma_records import provenance, freeze


def uniform_source_vacuum(radius, length, beam_radius):
    """Direct beam/direction quadrature, analytically averaged over initial z."""
    x, w = roots_legendre(96)
    rho = beam_radius*np.sqrt((x+1)/2)[:,None,None]
    mu = ((x+1)/2)[None,:,None]; phi = (np.pi*(x+1))[None,None,:]
    a = 1-mu*mu; b = rho*np.sqrt(a)*np.cos(phi)
    side = (-b+np.sqrt(b*b+a*(radius*radius-rho*rho)))/a
    zmax = length/mu
    side_probability = np.maximum(1-side/zmax, 0.)
    mean_flight = np.where(side <= zmax, side-side**2/(2*zmax), zmax/2)
    weights = w[:,None,None]*w[None,:,None]*w[None,None,:]/8
    return float(np.sum(weights*side_probability)), float(np.sum(weights*mean_flight)/c)


def run(destination, packets, batches):
    out = Path(destination)
    if out.exists() or packets < 8 or batches < 4:
        raise ValueError('Fresh output, >=8 histories and >=4 batches required')
    out.mkdir(parents=True); start = time.monotonic()
    identity = provenance(['examples.resonance_axial_experiment']); freeze(out, identity)
    root = Path(__file__).resolve().parents[1]
    shutil.copytree(root/'data/atomic/radiation_reference', out/'experimental_reference')
    shutil.copytree(root/'data/atomic/excitation', out/'atomic_inputs')
    model = RadiativeLevels.from_file(out/'atomic_inputs/argon_low_radiative.json')
    radius, length, beam, aperture, slit, temperature = .0111, .07, .0015, .00357, 50e-6, 293.
    references = {float(r['pressure_torr']):r for r in csv.DictReader(
        (out/'experimental_reference/table1.csv').open()) if r['later_impurity_flag']=='false'}
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
        atomic_input_sha256=model.input_sha256, histories_per_batch=packets, batches=batches,
        geometry=dict(radius_m=radius,length_m=length,beam_radius_m=beam,aperture_radius_m=aperture,slit_width_m=slit),
        temperature_k=temperature, source='Impulsive uniform axial source; uniform3mm beam disc assumed from collimator. Initial z integrated exactly with finite end losses.',
        detector='Parallel+x narrow view; actual conditional observer frequency, expected photons per steradian. No forced escape or amplitude fitting of simulated data.',
        collision_model='nu_c=0.904*A*(gu/gl)*N*lambda^3/(6*pi^2), ORNL-TM-4368 p35. Lorentz HWHM adds nu_c/(4*pi); reset probability conditional on Exp(A) hold t is1-exp(-nu_c*t). Full memory loss, not quenching.',
        limitations='Uniform stationary bath; scalar isotropic scattering, constant hnu0, no recoil/atom motion/quenching/optical spectral response. Actual beam and pulse profiles unknown. Three declared fit windows; published exact fit intervals unavailable. No coupled plasma or FRC formation claim.',
        completed_batches=[], comparisons=[])
    def save():
        report['elapsed_seconds'] = time.monotonic()-start
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    vacuum_edges = np.linspace(0, 30e-9, 151)
    vacuum = axial_transport(model, 0, 0., KB*temperature, radius, length, beam,
        aperture, slit, vacuum_edges, 65536, 20001)
    np.savez_compressed(out/'vacuum.npz', **vacuum)
    reference_hist, reference_integral = vacuum_reference(model, radius, beam, slit, length, vacuum_edges)
    side_reference, flight_reference = uniform_source_vacuum(radius, length, beam)
    exits = vacuum['exit_intervals']; weights = (exits[:,2]-exits[:,1])/length
    side_per_history = np.bincount(exits[:,0].astype(int), weights=weights*exits[:,5], minlength=65536)
    mean_exit = vacuum['mean_exit_time_per_history_s']
    report['vacuum'] = dict(histories=65536,
        detector_integral_per_sr=float(vacuum['detector_packet_integral_per_sr'].mean()),
        exact_detector_integral_per_sr=reference_integral,
        detector_cdf_maximum_difference=float(abs(np.cumsum(vacuum['detector_probability_per_sr']-reference_hist)/reference_integral).max()),
        side_fraction=float(side_per_history.mean()), side_fraction_standard_error=float(side_per_history.std(ddof=1)/256),
        reference_side_fraction=side_reference,
        mean_escape_ns=float(mean_exit.mean()*1e9), mean_escape_standard_error_ns=float(mean_exit.std(ddof=1)*1e9/256),
        reference_mean_escape_ns=float((1/model.A[0]+flight_reference)*1e9),
        maximum_source_weight_imbalance=float(abs(vacuum['escaped_source_weight']-1).max()))
    np.savez_compressed(out/'vacuum_analytic.npz', edges_s=vacuum_edges, probability_per_sr=reference_hist)
    save(); print('Integrated-source vacuum',report['vacuum'],flush=True)
    edges = np.linspace(0, 10e-6, 201)
    fig, axes = plt.subplots(2, 2, figsize=(12,8), constrained_layout=True)
    for pressure_index, pressure in enumerate((.0147, .042)):
        reference = references[pressure]; measured = float(reference['early_escape_rate_per_us'])
        density = pressure*(101325/760)/(KB*temperature)
        for mode_index, collisions in enumerate((False, True)):
            mode = 'resonant_collisions' if collisions else 'no_collisions'
            collision_rate = resonant_decorrelation_rate(model, 0, density) if collisions else 0.
            histograms=[]
            for batch in range(batches):
                seed = 30000+1000*pressure_index+100*mode_index+batch
                name=f'p{pressure:g}_{mode}_batch{batch:02}'
                def progress(generation, remaining):
                    if generation % 10000 == 0:
                        print(name,generation,'remaining',remaining,flush=True)
                result = axial_transport(model,0,density,KB*temperature,radius,length,beam,
                    aperture,slit,edges,packets,seed,decorrelation_rate_s=collision_rate,progress=progress)
                np.savez_compressed(out/f'{name}.npz',**result)
                histograms.append(result['detector_probability_per_sr'])
                row=dict(name=name,pressure_torr=pressure,mode=mode,seed=seed,
                    maximum_source_weight_imbalance=float(abs(result['escaped_source_weight']-1).max()),
                    mean_escape_us=float(result['mean_exit_time_per_history_s'].mean()*1e6),
                    detector_integral_per_sr=float(result['detector_packet_integral_per_sr'].mean()),
                    detector_outside_bins_per_sr=float(result['detector_outside_time_probability_per_sr']),
                    collision_rate_s=collision_rate,expected_collision_probability=float(result['collision_probability']),
                    observed_emission_collision_fraction=float(result['collision_emission_source_weight'].sum()/result['emission_source_weight'].sum()),
                    maximum_relative_detuning=float(result['maximum_relative_detuning']),
                    maximum_relative_observer_detuning=float(result['maximum_relative_observer_detuning']))
                report['completed_batches'].append(row);save()
                print('Complete',name,'detector/sr',row['detector_integral_per_sr'],flush=True)
            fits=summarize_batches(histograms,edges*1e6,measured)
            comparison=dict(pressure_torr=pressure,mode=mode,measured=reference,fits=fits)
            report['comparisons'].append(comparison);save();print('Comparison',comparison,flush=True)
            h=np.array(histograms);mean=h.mean(axis=0);sem=h.std(axis=0,ddof=1)/np.sqrt(batches)
            centres=.5*(edges[1:]+edges[:-1])*1e6;width=np.diff(edges)*1e6
            color=f'C{mode_index}';label='Resonant collisions' if collisions else 'No collisions'
            axes[0,pressure_index].errorbar(centres,mean/width,yerr=sem/width,fmt='.',markersize=2.5,
                                           color=color,alpha=.7,label=label+' ± batch SEM')
            if collisions:
                anchor=(centres>=.5)&(centres<3.);shape=np.exp(-measured*centres)
                amplitude=mean[anchor].sum()/np.sum(shape[anchor]*width[anchor])
                axes[0,pressure_index].plot(centres,amplitude*shape,'k--',label=f'Published slope {measured:g}/µs')
            axes[1,pressure_index].errorbar(np.arange(3)+(mode_index-.5)*.12,
                [f['decay_per_us'] for f in fits],yerr=[f['jackknife_standard_error_per_us'] for f in fits],
                fmt='o',color=color,label=label)
        ax=axes[0,pressure_index];ax.set(title=f'104.8 nm, {pressure*1000:g} mTorr',xlabel='Time [µs]',
            ylabel='Photons / initial quantum / sr / µs',yscale='log',xlim=(0,8));ax.legend(fontsize=8)
        ax=axes[1,pressure_index];error=measured*float(reference['reported_relative_error_bound'])
        ax.axhspan(measured-error,measured+error,color='black',alpha=.1,label='Published ± stated3% bound')
        ax.axhline(measured,color='black',linestyle='--')
        ax.set(xticks=np.arange(3),xticklabels=['0.25–2 µs','0.5–3 µs','1–5 µs'],
               xlabel='Fixed fit window',ylabel='Decay rate [1/µs]');ax.legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Argon radiation measured through a narrow side view\nFinite ends retained; uniform axial source integrated exactly; collision rates are sourced, not fitted')
    fig.savefig(out/'resonance_axial_experiment.png',dpi=170);fig.savefig(out/'resonance_axial_experiment.pdf');plt.close(fig)
    report['terminal']=True;save();print('Axial-source experiment complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination')
    parser.add_argument('--packets',type=int,default=4096);parser.add_argument('--batches',type=int,default=8)
    args=parser.parse_args();run(args.destination,args.packets,args.batches)
