"""Complete resonance-photon histories and independent free-flight comparisons."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, Boltzmann as KB
from scipy.integrate import quad
from scipy.special import roots_legendre, voigt_profile
from plasma_bench.line_radiation import RadiativeLevels
from plasma_bench.resonance_transport import transport
from .implicit_plasma_records import provenance, freeze


def centre_rays(radius, length):
    """Independent angular quadrature, split at the side/end transition."""
    half = length/2; cut = half/np.hypot(radius, half)
    nodes, weights = roots_legendre(96)
    mu = np.r_[cut*(nodes+1)/2, cut+(1-cut)*(nodes+1)/2]
    w = np.r_[weights*cut/2, weights*(1-cut)/2]
    distance = np.minimum(radius/np.sqrt(1-mu*mu), half/mu)
    return distance, w, cut


def first_escape_reference(model, line, density, temperature, radius, length):
    """Frequency/angle integral of Beer-Lambert transmission from the centre."""
    sigma, gamma = model.widths(temperature)
    damping = gamma[line]/sigma[line]
    strength = (density*model.A[line]*model.wavelength[line]**2/(8*np.pi)
        *model.weight[model.upper[line]]/model.weight[0]/sigma[line])
    distances, weights, _ = centre_rays(radius, length)
    def integrand(x):
        profile = voigt_profile(x, 1., damping)
        return 2*profile*(weights@np.exp(-strength*profile*distances))
    value, error = quad(integrand, 0, np.inf, epsabs=1e-10, epsrel=1e-9, limit=250)
    return value, error


def run(destination, packets):
    out = Path(destination)
    if out.exists():
        raise ValueError('A new output directory is required')
    out.mkdir(parents=True); start = time.monotonic()
    identity = provenance(['examples.resonance_transport']); freeze(out, identity)
    data = Path(__file__).resolve().parents[1]/'data/atomic/excitation'
    shutil.copytree(data, out/'atomic_inputs')
    model = RadiativeLevels.from_file(data/'argon_low_radiative.json')
    radius, length, temperature = .075, 1., KB*300.
    times = np.unique(np.r_[np.linspace(0, 100e-9, 201), np.linspace(0, 3.5e-6, 351),
                            np.geomspace(3.5e-6, 1e-3, 151)])
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
        atomic_input_sha256=model.input_sha256, radius_m=radius, length_m=length,
        atom_temperature_k=300., initial_source='Weak impulsive excitation at cylinder centre',
        scope='Uniform fixed neutral bath, scalar isotropic resonance transfer; coherent atom-frame redistribution compared with complete randomization at unchanged Voigt widths. Constant line-centre packet energy, exponential natural holding time. No collisions, ionization, recoil, atom displacement, polarization or surface return.',
        cases=[])
    def save():
        report['elapsed_seconds'] = time.monotonic()-start
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save(); fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    def case(line, pressure, count, redistribution, seed):
        density = pressure*(101325/760)*1e-3/(KB*300.)
        name = f'line{line}_p{pressure:g}_{redistribution}_n{count}'
        def progress(generation, remaining):
            print(name, 'generation', generation, 'remaining', remaining, flush=True)
        result = transport(model, line, density, temperature, radius, length, count,
                           seed, times, redistribution=redistribution, progress=progress)
        np.savez_compressed(out/f'{name}.npz', **result)
        reference, quadrature_error = first_escape_reference(model, line, density,
                                                             temperature, radius, length)
        escaped = result['exit_time_s']; unabsorbed = result['absorptions'] == 0
        mean = float(np.mean(escaped)); sem = float(np.std(escaped, ddof=1)/np.sqrt(count))
        entry = dict(name=name, line_index=line, wavelength_nm=float(model.wavelength[line]*1e9),
            neutral_pressure_equivalent_mtorr=pressure, neutral_density_m3=density,
            packets=count, seed=seed, redistribution=redistribution,
            mean_escape_time_ns=mean*1e9, mean_escape_standard_error_ns=sem*1e9,
            median_escape_time_ns=float(np.median(escaped)*1e9),
            quantile95_escape_time_ns=float(np.quantile(escaped,.95)*1e9),
            maximum_escape_time_ns=float(escaped.max()*1e9),
            mean_absorptions=float(result['absorptions'].mean()),
            maximum_absorptions=int(result['absorptions'].max()),
            unabsorbed_escape_fraction=float(unabsorbed.mean()),
            reference_unabsorbed_escape_fraction=reference,
            reference_frequency_quadrature_error=quadrature_error,
            unabsorbed_binomial_standard_error=float(np.sqrt(reference*(1-reference)/count)),
            escaped_by_3point5us_fraction=float(np.mean(escaped <= 3.5e-6)),
            mean_flight_time_ns=float(result['flight_time_s'].mean()*1e9),
            side_escape_fraction=float(result['exit_side'].mean()),
            maximum_packet_count_imbalance=int(abs(result['packet_balance']).max()),
            minimum_interior_count=int(min(result['excited_count'].min(),result['photon_count'].min())),
            maximum_relative_detuning=float(result['maximum_relative_detuning']))
        if pressure == 0:
            distances, weights, cut = centre_rays(radius, length)
            exact_cdf = np.array([weights@(-np.expm1(-model.A[line]*np.maximum(t-distances/c, 0))) for t in times])
            exact_mean = 1/model.A[line]+(radius*np.arcsin(cut)-length/2*np.log(cut))/c
            entry['vacuum_mean_escape_ns'] = exact_mean*1e9
            entry['vacuum_cdf_maximum_difference'] = float(abs(result['escaped_count']/count-exact_cdf).max())
            entry['vacuum_exact_side_fraction'] = cut
            np.savez_compressed(out/f'{name}_analytic.npz', time_s=times, escaped_fraction=exact_cdf)
        report['cases'].append(entry); save()
        print('Complete', name, 'mean ns', entry['mean_escape_time_ns'],
              '+/-', entry['mean_escape_standard_error_ns'],
              'escape by 3.5 us', entry['escaped_by_3point5us_fraction'], flush=True)
        if line == 1 and count == packets:
            label = f'{pressure:g} mTorr'
            if redistribution == 'atom_frame':
                axes[0,0].semilogx(times[1:]*1e9, 1-result['escaped_count'][1:]/count, label=label)
            if pressure in (1., 100.):
                axes[0,1].semilogx(times[1:]*1e9, 1-result['escaped_count'][1:]/count,
                    label=f'{label}, {redistribution.replace("_", " ")}')
        return entry
    for line in (0, 1):
        case(line, 0., packets, 'atom_frame', 1000+line)
        for index, pressure in enumerate((.1, 1., 10., 100.)):
            for mode_index, mode in enumerate(('atom_frame', 'complete')):
                case(line, pressure, packets, mode, 2000+100*line+10*index+mode_index)
    for count, seed in ((packets//4, 4101), (packets*4, 4102)):
        case(1, 10., count, 'atom_frame', seed)
    for line in (0, 1):
        for mode, style in (('atom_frame', '-'), ('complete', '--')):
            selected = [x for x in report['cases'] if x['line_index']==line and
                x['packets']==packets and x['neutral_pressure_equivalent_mtorr']>0 and x['redistribution']==mode]
            axes[1,0].errorbar([x['neutral_pressure_equivalent_mtorr'] for x in selected],
                [x['mean_escape_time_ns'] for x in selected],
                yerr=[2*x['mean_escape_standard_error_ns'] for x in selected],
                linestyle=style, marker='o', label=f'{model.wavelength[line]*1e9:.1f} nm {mode.replace("_", " ")}')
    for mode, marker in (('atom_frame','o'), ('complete','x')):
        selected=[x for x in report['cases'] if x['redistribution']==mode and x['neutral_pressure_equivalent_mtorr']>0]
        axes[1,1].errorbar([x['reference_unabsorbed_escape_fraction'] for x in selected],
            [x['unabsorbed_escape_fraction']-x['reference_unabsorbed_escape_fraction'] for x in selected],
            yerr=[2*x['unabsorbed_binomial_standard_error'] for x in selected],
            fmt=marker, label=mode.replace('_',' '))
    axes[0,0].set(title='106.7 nm energy still inside the gas', xlabel='Time [ns]', ylabel='Retained fraction', ylim=(0,1), xlim=(1,1e6))
    axes[0,1].set(title='Frequency correlation changes escape', xlabel='Time [ns]', ylabel='Retained fraction', ylim=(0,1), xlim=(1,1e6))
    axes[1,0].set(title='Mean escape time, bars = 2 standard errors', xlabel='Neutral pressure equivalent at 300 K [mTorr]', ylabel='Time [ns]', xscale='log', yscale='log')
    axes[1,1].axhline(0, color='grey', linewidth=1)
    axes[1,1].set(title='Independent first-flight transmission integral', xlabel='Exact unabsorbed escape probability', ylabel='Monte Carlo minus quadrature')
    for ax in axes.flat:
        ax.legend(fontsize=7); ax.grid(alpha=.2)
    fig.suptitle('Resonance photons move, excite atoms and eventually reach a boundary\nWeak centre pulse, uniform neutral argon, 1 m × 15 cm cylinder; transport component, not a formation prediction')
    fig.savefig(out/'resonance_transport.png', dpi=170); fig.savefig(out/'resonance_transport.pdf'); plt.close(fig)
    report['terminal'] = True; save(); print('Transport study complete', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--packets', type=int, default=4096)
    args = parser.parse_args(); run(args.destination, args.packets)
