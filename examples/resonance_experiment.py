"""Finite argon-cell photon histories and a published narrow-view decay signal."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann as KB, c
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri, roots_legendre
from plasma_bench.line_radiation import RadiativeLevels
from plasma_bench.resonance_transport import transport
from plasma_bench.resonance_detector import SideView
from .implicit_plasma_records import provenance, freeze


def beam_source(count, seed, radius, length):
    """Uniform physical beam disc/axis, importance sampled near the side view."""
    rng = np.random.default_rng(seed)
    radial = radius*np.sqrt(rng.random(count)); phi = rng.uniform(0, 2*np.pi, count)
    sigma = .0037; fraction_uniform = .25
    lo, hi = ndtr(-length/(2*sigma)), ndtr(length/(2*sigma))
    z = rng.uniform(-length/2, length/2, count)
    central = rng.random(count) >= fraction_uniform
    z[central] = sigma*ndtri(rng.uniform(lo, hi, np.sum(central)))
    proposal = (fraction_uniform/length+(1-fraction_uniform)
                *np.exp(-.5*(z/sigma)**2)/(np.sqrt(2*np.pi)*sigma*(hi-lo)))
    weights = 1/(length*proposal)
    return np.column_stack((radial*np.cos(phi), radial*np.sin(phi), z)), weights


def decay_rate(edges_us, histogram, window):
    """Fit exponential bin probabilities with a free overall count amplitude."""
    centres = .5*(edges_us[1:]+edges_us[:-1])
    use = (edges_us[:-1] >= window[0]-1e-10) & (edges_us[1:] <= window[1]+1e-10)
    t = centres[use]; h = histogram[use]
    if h.sum() <= 0 or len(t) < 3 or not np.allclose(np.diff(edges_us), np.diff(edges_us)[0]):
        raise ValueError('A populated interval and equal-width bins are required')
    measured = h@t/h.sum()
    def difference(rate):
        logp = -rate*t; p = np.exp(logp-logp.max())
        return p@t/p.sum()-measured
    # Negative apparent rates remain possible for an interval containing a rise.
    return float(brentq(difference, -100., 100.))


def summarize_batches(entries, edges_us, measured_rate):
    hist = np.asarray(entries); mean = hist.mean(axis=0)
    result = []
    for window in ((.25, 2.), (.5, 3.), (1., 5.)):
        fitted = decay_rate(edges_us, mean, window)
        jackknife = np.array([decay_rate(edges_us, np.delete(hist, k, axis=0).mean(axis=0), window)
                              for k in range(len(hist))])
        error = float(np.sqrt((len(hist)-1)/len(hist)*np.sum((jackknife-jackknife.mean())**2)))
        result.append(dict(window_us=list(window), decay_per_us=fitted,
                           jackknife_standard_error_per_us=error,
                           relative_difference_from_measured=fitted/measured_rate-1))
    return result


def vacuum_reference(model, radius, beam_radius, slit, length, edges):
    """Independent disc quadrature for an impulsive axial source in vacuum."""
    x, w = roots_legendre(72)
    radial = beam_radius*np.sqrt((x+1)/2)
    angle = np.pi*(x+1)
    xx = radial[:, None]*np.cos(angle)[None, :]
    yy = radial[:, None]*np.sin(angle)[None, :]
    flight = (np.sqrt(radius**2-yy**2)-xx)/c
    weights = w[:, None]*w[None, :]/4
    cdf = np.array([np.sum(weights*(-np.expm1(-model.A[0]*np.maximum(t-flight, 0)))) for t in edges])
    integral = slit/length/(4*np.pi)
    return np.diff(cdf)*integral, integral


def run(destination, packets, batches):
    out = Path(destination)
    if out.exists() or packets < 8 or batches < 4:
        raise ValueError('Fresh output, >=8 packets and >=4 independent batches required')
    out.mkdir(parents=True)
    start = time.monotonic(); identity = provenance(['examples.resonance_experiment'])
    freeze(out, identity)
    root = Path(__file__).resolve().parents[1]
    shutil.copytree(root/'data/atomic/radiation_reference', out/'experimental_reference')
    shutil.copytree(root/'data/atomic/excitation', out/'atomic_inputs')
    model = RadiativeLevels.from_file(out/'atomic_inputs/argon_low_radiative.json')
    radius, length, slit, aperture, beam_radius = .0111, .07, 50e-6, .00357, .0015
    rows = list(csv.DictReader((out/'experimental_reference/table1.csv').open()))
    references = {float(r['pressure_torr']):r for r in rows if r['later_impurity_flag']=='false'}
    edges = np.linspace(0, 10e-6, 201)
    times = np.unique(np.r_[edges, np.geomspace(10e-6, 1e-3, 101)])
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
        atomic_input_sha256=model.input_sha256, radius_m=radius, length_m=length,
        slit_width_m=slit, side_aperture_radius_m=aperture, beam_radius_m=beam_radius,
        temperature_k=300., physical_source='Uniform axial density and uniform beam disc, impulsive excitation; beam radius from collimator, actual profile not measured here.',
        source_sampling='q(z)=0.25/L+0.75 truncated Gaussian(sigma=3.7mm); true likelihood ratio (1/L)/q, never renormalized. Radial disc sampled directly.',
        detector='Parallel +x, constant small solid angle, circular aperture intersected by axial slit. Expected photons per steradian per initial physical quantum; no fitted optical amplitude in rate comparison.',
        limitations='Natural/thermal scalar atom-frame coherent model only. Fixed300K bath, no collisions/quenching/pressure broadening, atom motion, polarization, recoil or optical spectral response. Finite pulse width and actual beam profile unresolved. Published data fit windows are not numerically specified; three declared windows are retained. Not a formation calculation or demonstrated experimental agreement.',
        packets_per_batch=packets, batches=batches, completed_batches=[], comparisons=[])
    def save():
        report['elapsed_seconds'] = time.monotonic()-start
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    # An end-to-end vacuum trajectory independently checks view/arrival geometry.
    count = 262144
    positions, weights = beam_source(count, 9000, beam_radius, length)
    vacuum_edges = np.linspace(0, 30e-9, 151)
    observer = SideView(radius, aperture, slit, vacuum_edges, weights)
    result = transport(model, 0, 0., KB*300, radius, length, count, 9001,
                       vacuum_edges, initial_positions_m=positions, side_view=observer)
    result.update(observer.result()); np.savez_compressed(out/'vacuum.npz', **result)
    exact_hist, exact_integral = vacuum_reference(model, radius, beam_radius, slit, length, vacuum_edges)
    integral = observer.integral.mean(); error = observer.integral.std(ddof=1)/np.sqrt(count)
    report['vacuum'] = dict(packets=count, detector_integral_per_sr=float(integral),
        standard_error_per_sr=float(error), exact_integral_per_sr=exact_integral,
        standardized_integral_difference=float((integral-exact_integral)/error),
        maximum_packet_count_imbalance=int(abs(result['packet_balance']).max()))
    np.savez_compressed(out/'vacuum_analytic.npz', edges_s=vacuum_edges, probability_per_sr=exact_hist)
    save(); print('Vacuum view complete', report['vacuum'], flush=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for index, pressure in enumerate((.0147, .042)):
        reference = references[pressure]; measured = float(reference['early_escape_rate_per_us'])
        histograms = []
        for batch in range(batches):
            seed = 10000+1000*index+10*batch
            positions, weights = beam_source(packets, seed, beam_radius, length)
            observer = SideView(radius, aperture, slit, edges, weights)
            density = pressure*(101325/760)/(KB*300)
            name = f'p{pressure:g}_batch{batch:02}'
            def progress(generation, remaining):
                if generation % 10000 == 0:
                    print(name, generation, 'absorptions; remaining', remaining, flush=True)
            result = transport(model, 0, density, KB*300, radius, length, packets,
                seed+1, times, initial_positions_m=positions, side_view=observer, progress=progress)
            result.update(observer.result()); np.savez_compressed(out/f'{name}.npz', **result)
            histograms.append(result['detector_probability_per_sr'])
            entry = dict(name=name, pressure_torr=pressure, source_seed=seed, photon_seed=seed+1,
                mean_source_likelihood_ratio=float(weights.mean()),
                detector_integral_per_sr=float(observer.integral.mean()),
                detector_outside_bins_per_sr=float(observer.outside_time_bins/packets),
                maximum_packet_count_imbalance=int(abs(result['packet_balance']).max()),
                maximum_absorptions=int(result['absorptions'].max()),
                maximum_escape_time_us=float(result['exit_time_s'].max()*1e6),
                maximum_relative_detuning=float(result['maximum_relative_detuning']))
            report['completed_batches'].append(entry); save()
            print('Complete', name, 'detector probability/sr', entry['detector_integral_per_sr'], flush=True)
        fits = summarize_batches(histograms, edges*1e6, measured)
        comparison = dict(pressure_torr=pressure, measured=reference, fits=fits)
        report['comparisons'].append(comparison); save(); print('Comparison', comparison, flush=True)
        h = np.array(histograms); mean = h.mean(axis=0); sem = h.std(axis=0, ddof=1)/np.sqrt(batches)
        centres = .5*(edges[1:]+edges[:-1])*1e6; width = np.diff(edges)*1e6
        axes[index].errorbar(centres, mean/width, yerr=sem/width, fmt='.', markersize=3, alpha=.65,
                             label='Simulated view ± batch SEM')
        # One amplitude anchors the published SLOPE for visualization only.
        anchor = (centres >= .5) & (centres < 3.)
        shape = np.exp(-measured*centres)
        amplitude = mean[anchor].sum()/np.sum(shape[anchor]*width[anchor])
        axes[index].plot(centres, amplitude*shape, label=f'Published slope {measured:g}/µs')
        axes[index].set(title=f'104.8 nm, {pressure*1000:g} mTorr', xlabel='Time [µs]',
                        ylabel='Photons / initial quantum / sr / µs', yscale='log', xlim=(0,8))
        axes[index].legend(fontsize=8); axes[index].grid(alpha=.2)
        for j, fit in enumerate(fits):
            axes[2].errorbar(index+(j-1)*.15, fit['decay_per_us'],
                yerr=fit['jackknife_standard_error_per_us'], fmt='o', color=f'C{j}',
                label=f"{fit['window_us'][0]:g}–{fit['window_us'][1]:g} µs fit" if index==0 else None)
        axes[2].errorbar(index, measured, yerr=measured*float(reference['reported_relative_error_bound']),
                         fmt='s', color='black', label='Published ± stated bound' if index==0 else None)
    axes[2].set(title='Decay rate depends on observation window', ylabel='Decay rate [1/µs]',
                xticks=[0,1], xticklabels=['14.7 mTorr','42 mTorr'])
    axes[2].legend(fontsize=8); axes[2].grid(alpha=.2)
    fig.suptitle('Finite argon cell, axial beam and narrow side view\nNatural/thermal correlated scattering; omitted broadening/quenching and source/optical uncertainties remain')
    fig.savefig(out/'resonance_experiment.png', dpi=170); fig.savefig(out/'resonance_experiment.pdf'); plt.close(fig)
    report['terminal'] = True; save(); print('Experimental view study complete', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--packets', type=int, default=4096)
    parser.add_argument('--batches', type=int, default=8)
    args = parser.parse_args(); run(args.destination, args.packets, args.batches)
