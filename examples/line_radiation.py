"""Complete resolved argon cascades and absorption at provisional tube scales."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV, Boltzmann as KB
from plasma_bench.line_radiation import RadiativeLevels
from .implicit_plasma_records import provenance, freeze


def evolve(model, initial, duration, steps):
    """Integrate actual levels and each line's cumulative photon count with RK4."""
    count = len(model.energy)
    y = np.r_[initial, np.zeros(len(model.A))]
    times = np.linspace(0, duration, steps+1)
    states = [y.copy()]
    def rhs(state):
        df, photons = model.source(state[:count])
        return np.r_[df, photons]
    for begin, end in zip(times[:-1], times[1:]):
        dt = end-begin
        k1 = rhs(y); k2 = rhs(y+dt*k1/2)
        k3 = rhs(y+dt*k2/2); k4 = rhs(y+dt*k3)
        y += dt*(k1+2*k2+2*k3+k4)/6
        model.source(y[:count])
        states.append(y.copy())
    states = np.array(states)
    stored = states[:, :count]@model.energy
    emitted = states[:, count:]@model.photon_energy
    return dict(time_s=times, populations=states[:, :count],
                photons_per_line=states[:, count:], stored_energy_j=stored,
                emitted_energy_j=emitted, energy_balance_j=stored+emitted-stored[0])


def sequential_reference(times):
    """Bateman solution for actual Ar 2p5 -> 1s4 -> ground decay."""
    # Independent scalar solution, directly using the published probabilities.
    a, b = 4.0e7, 1.320e8
    upper = np.exp(-a*times)
    middle = a/(b-a)*(np.exp(-a*times)-np.exp(-b*times))
    lower = 1-upper-middle
    return np.array([lower, middle, upper]).T


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('A new destination is required')
    out.mkdir(parents=True); started = time.monotonic()
    identity = provenance(['examples.line_radiation']); freeze(out, identity)
    data = Path(__file__).resolve().parents[1]/'data/atomic/excitation'
    shutil.copytree(data, out/'atomic_inputs')
    model = RadiativeLevels.from_file(data/'argon_low_radiative.json')
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
                  atomic_input_sha256=model.input_sha256, reference=[], cascade=[],
                  scope='Finite low-level spontaneous cascades only. No collisions or transfer. Produced photons equal escaped photons only in the explicit optically thin limit; opacity calculation separately evaluates that limit.')
    def save():
        report['elapsed_seconds'] = time.monotonic()-started
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for steps in (50, 100, 200):
        initial = np.zeros(15); initial[10] = 1.
        result = evolve(model, initial, 50e-9, steps)
        exact = sequential_reference(result['time_s'])
        error = abs(result['populations'][:, [0, 2, 10]]-exact)
        result['analytic_populations_ground_1s4_2p5'] = exact
        np.savez_compressed(out/f'sequential_t{steps}.npz', **result)
        entry = dict(steps=steps, maximum_population_error=float(error.max()),
                     maximum_energy_error_ev=float(abs(result['energy_balance_j']).max()/EV))
        report['reference'].append(entry); print('Sequential cascade', entry, flush=True)
        axes[1, 0].semilogy(result['time_s'][1:]*1e9, error[1:].max(axis=1), label=f'{steps} steps')
    previous = None
    for steps in (3500, 7000, 14000):
        initial = np.zeros(15); initial[0] = .9; initial[5:] = .01
        result = evolve(model, initial, 3.5e-6, steps)
        np.savez_compressed(out/f'argon_cascade_t{steps}.npz', **result)
        entry = dict(steps=steps, dt_ns=3500/steps,
                     final_populations=result['populations'][-1].tolist(),
                     final_metastable_fraction=float(result['populations'][-1, [1, 3]].sum()),
                     initial_stored_ev=float(result['stored_energy_j'][0]/EV),
                     final_stored_ev=float(result['stored_energy_j'][-1]/EV),
                     final_produced_photon_ev=float(result['emitted_energy_j'][-1]/EV),
                     minimum_population=float(result['populations'].min()),
                     maximum_number_error=float(abs(result['populations'].sum(axis=1)-1).max()),
                     maximum_energy_error_ev=float(abs(result['energy_balance_j']).max()/EV))
        if previous is not None:
            entry['maximum_population_change_from_double_dt'] = float(abs(previous['populations']-result['populations'][::2]).max())
            entry['maximum_emitted_energy_change_from_double_dt_ev'] = float(abs(previous['emitted_energy_j']-result['emitted_energy_j'][::2]).max()/EV)
        report['cascade'].append(entry); previous = result; save()
        print('Full argon cascade', entry, flush=True)
    t = result['time_s']*1e9; pop = result['populations']
    axes[0, 0].plot(t, pop[:, 5:].sum(axis=1), label='4p states')
    axes[0, 0].plot(t, pop[:, [2, 4]].sum(axis=1), label='4s resonance states')
    axes[0, 0].plot(t, pop[:, [1, 3]].sum(axis=1), label='4s metastable states')
    axes[0, 0].set(xlim=(0, 350), title='Finite cascade retains metastable atoms', xlabel='Time [ns]', ylabel='Fraction of neutral atoms')
    axes[0, 1].plot(t, result['stored_energy_j']/EV, label='Stored excitation')
    axes[0, 1].plot(t, result['emitted_energy_j']/EV, label='Produced photons')
    axes[0, 1].set(xlim=(0, 350), title='Atomic energy is not all emitted', xlabel='Time [ns]', ylabel='Energy [eV / neutral atom]')
    axes[1, 0].set(title='Independent two-stage analytic solution', xlabel='Time [ns]', ylabel='Largest population error')

    # These are ground-state neutral columns, not an inferred ionization history.
    pressure_mtorr = np.geomspace(.1, 100., 121)
    neutral_density = pressure_mtorr*(101325/760)*1e-3/(KB*300.)
    ground = np.zeros(15); ground[0] = 1.
    coefficient_per_density = model.absorption(ground, KB*300.)
    resonance = np.flatnonzero(model.lower == 0)
    radius = .075
    tau = neutral_density[:, None]*radius*coefficient_per_density
    np.savez_compressed(out/'cold_ground_neutral_columns.npz',
        pressure_mtorr=pressure_mtorr, neutral_density_m3=neutral_density,
        line_centre_tau=tau, wavelength_m=model.wavelength, path_m=radius)
    report['opacity'] = dict(atom_temperature_k=300., pressure_temperature_k=300.,
        neutral_fraction=1., radial_path_m=radius, axial_half_length_m=.5,
        scope='Uniform stationary ground-state neutral gas, thermal plus natural widths. Path is centre to wall/end; not a volume-averaged escape probability or transfer solution.', resonance_lines=[])
    for line in resonance:
        absorption = coefficient_per_density[line]
        wavelength = model.wavelength[line]*1e9
        entry = dict(line_index=int(line), vacuum_wavelength_nm=float(wavelength),
            absorption_cross_section_line_centre_m2=float(absorption),
            tau_radial_at_point1_mtorr=float(tau[0, line]),
            tau_radial_at_100_mtorr=float(tau[-1, line]),
            ground_neutral_density_for_tau1_m3=float(1/(radius*absorption)),
            natural_lifetime_ns=float(1e9/model.decay_rate[model.upper[line]]))
        report['opacity']['resonance_lines'].append(entry)
        axes[1, 1].loglog(pressure_mtorr, tau[:, line], label=f'{wavelength:.3f} nm')
    axes[1, 1].axhline(1, color='grey', linestyle='--', label='Optical depth 1')
    axes[1, 1].set(title='Resonance photons encounter neutral absorption', xlabel='Ground-neutral pressure equivalent at 300 K [mTorr]', ylabel='Line-centre optical depth over 7.5 cm')
    # Same nucleus density as the earlier homogeneous DI/RR source examples.
    source_density = 8.961429541213399e19
    report['source_density_neutral_columns'] = []
    for temperature in (KB*300., EV):
        for fraction in (.001, .01, .1, 1.):
            n = ground*source_density*fraction
            kappa = model.absorption(n, temperature)
            report['source_density_neutral_columns'].append(dict(
                nucleus_density_m3=source_density, ground_neutral_fraction=fraction,
                atom_temperature_ev=temperature/EV,
                radial_tau_resonance=(radius*kappa[resonance]).tolist(),
                axial_tau_resonance=(.5*kappa[resonance]).tolist()))
    for ax in axes.flat:
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Resolved neutral-argon radiative physics\nSpontaneous-only cascade and separate absorption calculation; no collisional-radiative or formation prediction')
    fig.savefig(out/'argon_line_radiation.png', dpi=170)
    fig.savefig(out/'argon_line_radiation.pdf'); plt.close(fig)
    report['terminal'] = True; save()
    print('Resolved radiation study complete', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    run(parser.parse_args().destination)
