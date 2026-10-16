"""Polarized finite-cell resonance histories compared with raw measured counts."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as KB
from scipy.optimize import brentq, minimize_scalar
from plasma_bench.line_radiation import RadiativeLevels
from plasma_bench.resonance_axial import axial_transport, resonant_decorrelation_rate
from plasma_bench.resonance_polarization import dipole_emission
from .resonance_axial_experiment import uniform_source_vacuum
from .implicit_plasma_records import provenance, freeze


def fit_intensity(counts, shape, background):
    """Profile only unknown source/detector intensity; keep the shape fixed."""
    counts, shape = np.asarray(counts), np.asarray(shape)
    if background == 0:
        return float(counts.sum()/shape.sum())
    def score(amplitude):
        return np.sum(shape*(counts/(background+amplitude*shape)-1))
    if score(0.) <= 0:
        return 0.
    return float(brentq(score, 0., counts.sum()/shape.sum(), xtol=1e-14))


def fit_exponential(counts, width_us, background):
    """Poisson exponential fit; equal-width integration is absorbed in amplitude."""
    counts = np.asarray(counts, float); times = np.arange(len(counts))*width_us
    def objective(log_rate):
        shape = np.exp(-np.exp(log_rate)*times)
        amplitude = fit_intensity(counts, shape, background)
        mu = background+amplitude*shape
        if np.any(mu <= 0):
            return np.inf
        return float(np.sum(mu-counts*np.log(mu))/counts.sum())
    solution = minimize_scalar(objective, bounds=(-9., 5.), method='bounded', options={'xatol':1e-10})
    if not solution.success or not -8.99 < solution.x < 4.99:
        raise RuntimeError('Decay fit failed or reached its numerical search limits')
    rate = float(np.exp(solution.x)); shape = np.exp(-rate*times)
    amplitude = fit_intensity(counts, shape, background); signal = amplitude*shape
    mu = background+signal; derivative = np.column_stack((signal, -rate*times*signal))
    covariance = np.linalg.inv(derivative.T@(derivative/mu[:,None]))
    return dict(rate_per_us=rate, conditional_standard_error_per_us=float(rate*np.sqrt(covariance[1,1])),
                signal_first_bin=amplitude)


def channel_sums(histograms, bins, offset_half_channels):
    start = offset_half_channels
    return np.asarray(histograms)[:, start:start+2*bins].reshape(-1,bins,2).sum(axis=2)


def run(destination, packets, batches):
    out = Path(destination)
    if out.exists() or packets < 8 or batches < 4:
        raise ValueError('Fresh output, >=8 histories per batch and >=4 batches required')
    out.mkdir(parents=True); started = time.monotonic()
    identity = provenance(['examples.resonance_channel_experiment']); freeze(out, identity)
    root = Path(__file__).resolve().parents[1]
    shutil.copytree(root/'data/atomic/radiation_reference',out/'experimental_reference')
    shutil.copytree(root/'data/atomic/excitation',out/'atomic_inputs')
    metadata = json.loads((out/'experimental_reference/channel_metadata.json').read_text())
    raw = list(csv.DictReader((out/'experimental_reference/channel_counts.csv').open()))
    model = RadiativeLevels.from_file(out/'atomic_inputs/argon_low_radiative.json')
    radius,length,beam,aperture,slit,temperature = .0111,.07,.0015,.00357,50e-6,293.
    report = dict(terminal=False, source_sha256=identity['source_sha256'], atomic_input_sha256=model.input_sha256,
        histories_per_batch=packets,batches=batches,temperature_k=temperature,
        geometry=dict(radius_m=radius,length_m=length,beam_radius_m=beam,aperture_radius_m=aperture,slit_width_m=slit),
        scope='Complete weak stationary-bath radiation histories; unaligned source, isolated zero-field dipoles, full orientation reset at collisions. No fitted collision rate, moving plasma coupling or claim of experimental validation.',
        inference='Fixed reported background, all three within-channel offsets. Amplitude-only count comparison; conditional Fisher data errors and whole-batch jackknife Monte Carlo errors. Conditional deviance is not a calibrated p-value.',
        completed_batches=[],comparisons=[])
    def save():
        report['elapsed_seconds']=time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    # Independent angular moments; no finite-cylinder implementation enters their references.
    rng=np.random.default_rng(41001);n=131072;dipoles=np.tile([1.,0.,0.],(n,1))
    directions,polarization=dipole_emission(rng,dipoles);moments=directions**2
    report['dipole_angular_limit']=dict(second_moments=moments.mean(axis=0).tolist(),
        standard_errors=(moments.std(axis=0,ddof=1)/np.sqrt(n)).tolist(),reference=[.2,.4,.4],
        maximum_transversality=float(abs(np.sum(directions*polarization,axis=1)).max()))
    phi=rng.uniform(0,2*np.pi,n);dipoles=np.column_stack((np.cos(phi),np.sin(phi),np.zeros(n)))
    directions,_=dipole_emission(rng,dipoles);mu2=directions[:,2]**2
    report['unpolarized_beam_limit']=dict(mean_mu_squared=float(mu2.mean()),standard_error=float(mu2.std(ddof=1)/np.sqrt(n)),reference=.4)
    vacuum=axial_transport(model,0,0.,KB*temperature,radius,length,beam,aperture,slit,
        np.linspace(0,30e-9,151),65536,41002,scattering='dipole')
    np.savez_compressed(out/'vacuum.npz',**vacuum)
    side_reference,flight_reference=uniform_source_vacuum(radius,length,beam)
    exits=vacuum['exit_intervals'];weights=(exits[:,2]-exits[:,1])/length
    side=np.bincount(exits[:,0].astype(int),weights=weights*exits[:,5],minlength=65536)
    scores=vacuum['detector_packet_stokes_integral_per_sr'];mean_exit=vacuum['mean_exit_time_per_history_s']
    report['vacuum']=dict(detector_stokes_mean_per_sr=scores.mean(axis=0).tolist(),
        detector_stokes_standard_error_per_sr=(scores.std(axis=0,ddof=1)/256).tolist(),
        reference_stokes_per_sr=[slit/(4*np.pi*length),0.,0.],side_fraction=float(side.mean()),
        side_standard_error=float(side.std(ddof=1)/256),reference_side_fraction=side_reference,
        mean_escape_ns=float(mean_exit.mean()*1e9),mean_escape_standard_error_ns=float(mean_exit.std(ddof=1)*1e9/256),
        reference_mean_escape_ns=float((1/model.A[0]+flight_reference)*1e9),
        maximum_source_weight_imbalance=float(abs(vacuum['escaped_source_weight']-1).max()))
    save();print('Complete dipole angular and vacuum limits',flush=True)
    fig,axes=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    for pindex,reference in enumerate(metadata['runs']):
        dt=reference['channel_width_us'];pressure=reference['pressure_torr'];period=reference['pulse_period_us']*1e-6
        edges=np.r_[np.arange(0,period,dt*.5e-6),period]
        counts={int(r['channel']):int(r['counts']) for r in raw if int(r['run'])==reference['run']}
        observed=np.array([counts[reference['first_full_channel']-j] for j in range(24)])
        background=reference['photon_background_per_s']*reference['acquisition_s']*dt/reference['pulse_period_us']
        windows=(6,9,15,24);data_fits=[dict(bins=nb,**fit_exponential(observed[:nb],dt,background)) for nb in windows]
        data=dict(run=reference['run'],pressure_torr=pressure,background_per_channel=background,
                  first_full_channel=reference['first_full_channel'],measured_fits=data_fits,models=[])
        report['comparisons'].append(data)
        axes[0,pindex].errorbar((np.arange(24)+1)*dt,observed-background,yerr=np.sqrt(observed),
            fmt='o',color='black',markersize=3,label='Measured counts; reported background subtracted')
        axes[1,pindex].errorbar(np.arange(4),[f['rate_per_us'] for f in data_fits],
            yerr=[f['conditional_standard_error_per_us'] for f in data_fits],fmt='ko',label='Raw-count Poisson fit')
        density=pressure*(101325/760)/(KB*temperature);nu=resonant_decorrelation_rate(model,0,density)
        for mindex,mode in enumerate(('isotropic','dipole')):
            histograms=[]
            for batch in range(batches):
                seed=42000+1000*pindex+100*mindex+batch;name=f'run{reference["run"]}_{mode}_{batch:02}'
                result=axial_transport(model,0,density,KB*temperature,radius,length,beam,aperture,slit,
                    edges,packets,seed,decorrelation_rate_s=nu,scattering=mode,
                    detector_period_s=period,spectral_thresholds_relative=(1e-4,1e-3,1e-2))
                np.savez_compressed(out/f'{name}.npz',**result);histograms.append(result['detector_probability_per_sr'])
                row=dict(name=name,seed=seed,mode=mode,pressure_torr=pressure,
                    maximum_source_weight_imbalance=float(abs(result['escaped_source_weight']-1).max()),
                    mean_escape_us=float(result['mean_exit_time_per_history_s'].mean()*1e6),
                    collision_rate_s=nu,detector_integral_per_sr=float(result['detector_packet_integral_per_sr'].mean()),
                    detector_stokes_integral_per_sr=result['detector_packet_stokes_integral_per_sr'].mean(axis=0).tolist(),
                    spectral_tail_integrals_per_sr=result['detector_spectral_tail_integral_per_sr'].tolist(),
                    maximum_relative_detuning=float(result['maximum_relative_detuning']),
                    maximum_relative_observer_detuning=float(result['maximum_relative_observer_detuning']))
                report['completed_batches'].append(row);save();print('Complete',name,flush=True)
            h=np.array(histograms);fits=[]
            for nb in windows:
                offsets=[]
                for offset in (0,1,2):
                    bins=channel_sums(h,nb,offset);mean=bins.mean(axis=0)
                    rate=fit_exponential(mean,dt,0.)['rate_per_us']
                    jack=np.array([fit_exponential((bins.sum(axis=0)-one)/(batches-1),dt,0.)['rate_per_us'] for one in bins])
                    scale=fit_intensity(observed[:nb],mean,background);prediction=background+scale*mean
                    deviance=2*np.sum(prediction-observed[:nb]+observed[:nb]*np.log(observed[:nb]/prediction))
                    offsets.append(dict(offset_us=offset*dt/2,rate_per_us=rate,
                        jackknife_standard_error_per_us=float(np.sqrt((batches-1)/batches*np.sum((jack-jack.mean())**2))),
                        amplitude_scale=scale,conditional_poisson_deviance=float(deviance),predicted_counts=prediction.tolist()))
                fits.append(dict(bins=nb,offsets=offsets))
            data['models'].append(dict(scattering=mode,fits=fits));save()
            bins=channel_sums(h,24,1);mean=bins.mean(axis=0);scale=fits[-1]['offsets'][1]['amplitude_scale']
            sem=bins.std(axis=0,ddof=1)/np.sqrt(batches);t=(np.arange(24)+1)*dt
            axes[0,pindex].plot(t,scale*mean,color=f'C{mindex}',label=mode+'; intensity fitted, rates fixed')
            axes[0,pindex].fill_between(t,scale*(mean-sem),scale*(mean+sem),color=f'C{mindex}',alpha=.18)
            axes[1,pindex].errorbar(np.arange(4)+(mindex+1)*.1,[f['offsets'][1]['rate_per_us'] for f in fits],
                yerr=[f['offsets'][1]['jackknife_standard_error_per_us'] for f in fits],fmt='o',color=f'C{mindex}',label=mode+'; midpoint offset')
            print('Compared',reference['run'],mode,flush=True)
        axes[0,pindex].set(title=f'Run {reference["run"]}, {pressure*1000:g} mTorr',xlabel='Time [µs], nominal half-channel offset',ylabel='Counts per channel above reported background')
        axes[1,pindex].set(xticks=np.arange(4),xticklabels=[f'{v*dt:.2f}' for v in windows],xlabel='Full fit-window length [µs]',ylabel='Decay rate [1/µs]')
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Argon resonance radiation against original detector counts\nSourced collisions; dipole polarization memory versus scalar control; time-zero uncertainty retained')
    fig.savefig(out/'channel_comparison.png',dpi=170);plt.close(fig)
    report['terminal']=True;save();print('Channel experiment complete',out,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination')
    parser.add_argument('--packets',type=int,default=4096);parser.add_argument('--batches',type=int,default=8)
    args=parser.parse_args();run(args.destination,args.packets,args.batches)
