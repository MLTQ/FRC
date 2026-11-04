"""Measure post-pulse ion response from complete native states, without evolution."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull


def read_json(path):
    return json.loads(path.read_text())


def first_crossing(times, values, target, start, *, below=False):
    """Return an observed sample bracket, never an interpolated physical state."""
    selected = np.flatnonzero((times >= start) & ((values <= target) if below else (values >= target)))
    if not len(selected):
        return dict(reached=False, observed_through_us=float(times[-1]),
                    post_release_window_us=float(times[-1]-start), target=float(target))
    i = int(selected[0])
    previous = max(i-1, int(np.searchsorted(times, start)))
    return dict(reached=True, target=float(target), time_us=float(times[i]),
                since_release_us=float(times[i]-start),
                bracket_us=[float(times[previous]), float(times[i])])


def compare_archives(run, reference, output):
    """Compare every native array at each shared original archive time."""
    dense, coarse = read_json(run/'history.json'), read_json(reference/'history.json')
    lookup = {round(row['time_ns'], 7): i for i, row in enumerate(dense)}
    comparisons = []
    for j, row in enumerate(coarse):
        i = lookup[round(row['time_ns'], 7)]
        with np.load(run/f'state_{i:04d}.npz') as a, np.load(reference/f'state_{j:04d}.npz') as b:
            same = set(a.files) == set(b.files)
            keys = {key: np.array_equal(a[key], b[key]) for key in b.files}
        comparisons.append(dict(time_ns=row['time_ns'], same_keys=same, arrays=keys))
    result = dict(shared_states=len(comparisons),
                  all_bitwise_equal=all(c['same_keys'] and all(c['arrays'].values()) for c in comparisons),
                  scope='Archive-cadence independence, not physical convergence', states=comparisons)
    (output/'shared_times.json').write_text(json.dumps(result, indent=2)+'\n')
    return {key: value for key, value in result.items() if key != 'states'}


def measure(run, output, reference=None):
    output.mkdir(parents=True, exist_ok=False)
    manifest, history = read_json(run/'run.json'), read_json(run/'history.json')
    geometry = read_json(run/'geometry.json')
    vertices = np.asarray(geometry['vertices'])
    # The closed hull is only an occupancy diagnostic; its open faces are not metal.
    equations = np.unique(np.round(ConvexHull(vertices).equations, 13), axis=0)
    outer_radius = float(np.linalg.norm(vertices, axis=1).max()+geometry['wire_radius'])
    switches = read_json(run/'drive_switches.json')
    release = next(event['time_us'] for event in switches if not event['on'])
    times = np.array([row['time_us'] for row in history])
    release_index = int(np.argmin(abs(times-release)))
    initial = history[0]['ions']
    rows, max_radius_error, max_speed_error = [], 0., 0.
    for i, saved in enumerate(history):
        with np.load(run/f'state_{i:04d}.npz') as state:
            select = state['species'] == 1
            x, v = state['x'][select], state['v'][select]
            if not len(x):
                raise ValueError('This response analysis requires a nonempty ion population')
            if not (np.isfinite(x).all() and np.isfinite(v).all()):
                raise ValueError(f'Nonfinite native ion state {i}')
        r = np.linalg.norm(x, axis=1)
        dx, dv = x-x.mean(axis=0), v-v.mean(axis=0)
        radius = float(np.sqrt(np.mean(np.sum(dx*dx, axis=1))))
        speed = float(np.mean(np.sum(dx*dv, axis=1))/radius)
        max_radius_error = max(max_radius_error, abs(radius-saved['ion_motion']['radius_rms_m']))
        max_speed_error = max(max_speed_error, abs(speed-saved['ion_motion']['radial_speed_m_s']))
        hull_inside = np.all(x@equations[:, :3].T+equations[:, 3] <= 1e-12, axis=1)
        absorbed = np.asarray(saved['collected_macroparticles'])[:, 1]
        if len(x)+int(absorbed.sum()) != initial:
            raise ValueError(f'Ion count ledger does not close at state {i}')
        q10, q50, q90 = np.quantile(r, [.1, .5, .9])
        rows.append(dict(time_us=saved['time_us'], since_release_us=saved['time_us']-release,
                         voltage_v=saved['voltage_v'], centre_potential_v=saved['centre_potential_v'],
                         ions=len(x), inside_hull=int(hull_inside.sum()),
                         inside_outer_sphere=int((r < outer_radius).sum()),
                         core_ions=int((r < .02).sum()), absorbed_cage=int(absorbed[0]),
                         absorbed_wall=int(absorbed[1]), radius_rms_m=radius,
                         origin_r10_m=float(q10), origin_r50_m=float(q50), origin_r90_m=float(q90),
                         radial_speed_m_s=speed,
                         outward_moving_fraction=float(np.mean(np.sum(dx*dv, axis=1) > 0)),
                         mean_energy_ev=saved['ion_motion']['mean_energy_ev']))
    with (output/'time_series.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    columns = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    post = times >= release
    at_release = rows[release_index]
    occupied = columns['inside_hull']
    baseline = at_release['inside_hull']
    response = {}
    for fraction in (.9, .8, .5, 1/np.e):
        response[f'cage_occupancy_fraction_{fraction:.6f}'] = first_crossing(
            times, occupied, baseline*fraction, release, below=True)
    for fraction in (1.1, 1.25, 1.5):
        response[f'radius_factor_{fraction}'] = first_crossing(
            times, columns['radius_rms_m'], at_release['radius_rms_m']*fraction, release)
    response['first_ion_beyond_outer_radius'] = first_crossing(
        times, columns['ions']-columns['inside_outer_sphere'], 1, release)
    response['core_population_half'] = first_crossing(
        times, columns['core_ions'], at_release['core_ions']*.5, release, below=True)
    response['voltage_below_one_percent'] = first_crossing(
        times, columns['voltage_v'], at_release['voltage_v']*.01, release, below=True)
    # Fit the prescribed first-order voltage fall; this is not an ion relaxation fit.
    fall = (times > release) & (times <= release+.2)
    slope = np.polyfit(times[fall]-release, np.log(columns['voltage_v'][fall]), 1)[0]
    always_outward = bool(np.all(columns['radial_speed_m_s'][post] > 0))
    conclusion = ('Instantaneous RMS dilation remains outward at every post-release sample; '
                  'no bounded breathing cycle is identified.' if always_outward else
                  'Dilation changes sign or vanishes; inspect phase space before inferring a cycle.')
    if not response['cage_occupancy_fraction_0.500000']['reached']:
        conclusion += ' Cage half-empty threshold is not reached within this record.'
    if not response['cage_occupancy_fraction_0.367879']['reached']:
        conclusion += ' Cage e-fold occupancy threshold is not reached within this record.'
    result = dict(run=str(run.resolve()), status=manifest['status'], states=len(rows),
                  release_us=release, observed_until_us=float(times[-1]),
                  post_release_window_us=float(times[-1]-release),
                  definitions=dict(cage='Convex hull of wire centreline vertices; open-face diagnostic only',
                                   outer_radius_m=outer_radius,
                                   radius='RMS about the instantaneous live-ion centroid; losses can bias it',
                                   crossing='First saved sample meeting threshold, with previous-sample bracket',
                                   survival='Live anywhere in chamber, distinct from cage occupancy'),
                  initial=rows[0], at_release=at_release, final=rows[-1], response=response,
                  voltage_fall_tau_ns=float(-1000/slope),
                  post_release_radial_speed_min_m_s=float(columns['radial_speed_m_s'][post].min()),
                  final_instantaneous_expansion_time_us=float(rows[-1]['radius_rms_m']/rows[-1]['radial_speed_m_s']*1e6),
                  relaxation_fit=None,
                  conclusion=conclusion,
                  limits='Collisionless development model; no Coulomb/neutral thermal relaxation, unresolved force/phase/mesh/particle and near-wire-loss convergence. Thresholds describe this realization, not a validated confinement prediction.',
                  native_checks=dict(max_radius_vs_saved_error_m=max_radius_error,
                                     max_speed_vs_saved_error_m_s=max_speed_error,
                                     ion_count_ledger_closes_all_states=True),
                  inputs_sha256={name: hashlib.sha256((run/name).read_bytes()).hexdigest()
                                 for name in ('run.json', 'history.json', 'geometry.json', 'drive_switches.json')})
    if reference:
        result['shared_time_comparison'] = compare_archives(run, reference, output)
    (output/'relaxation.json').write_text(json.dumps(result, indent=2)+'\n')
    plot(columns, release, outer_radius, initial, output)
    return result


def plot(c, release, outer_radius, initial, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.6), constrained_layout=True)
    t = c['time_us']
    a = axes[0, 0]
    a.plot(t, c['radius_rms_m']*100, label='RMS about ion centroid', color='#167f91', lw=2)
    a.plot(t, c['origin_r50_m']*100, label='Median distance from centre', color='#cf8124')
    a.plot(t, c['origin_r90_m']*100, label='90th percentile distance', color='#8368b1')
    a.axhline(outer_radius*100, color='#555555', ls=':', label='Cage outer radius')
    a.set(ylabel='Ion cloud radius (cm)', title='Cloud expands after the pulse')
    a.legend(fontsize=8)
    a = axes[0, 1]
    a.plot(t, c['inside_hull']/initial*100, label='Inside cage envelope', color='#167f91', lw=2)
    a.plot(t, c['ions']/initial*100, label='Live anywhere in chamber', color='#cf8124')
    a.set(ylabel='Fraction of initial ions (%)', ylim=(0, 104), title='Inside the cage ≠ still in the chamber')
    a.legend(fontsize=8, loc='lower left')
    a = axes[1, 0]
    a.plot(t, c['radial_speed_m_s']/1000, color='#167f91', lw=2)
    a.axhline(0, color='#555555', ls=':')
    a.set(ylabel='RMS expansion speed (km/s)', title='Mean expansion stays outward')
    a = axes[1, 1]
    a.plot(t, c['voltage_v'], label='Cage voltage', color='#cf8124')
    a.plot(t, c['centre_potential_v'], label='Centre potential', color='#167f91')
    a.set(yscale='symlog', ylim=(-.05, 5000), ylabel='Potential relative to wall (V)',
          title='Applied pulse falls; space charge remains')
    a.legend(fontsize=8)
    for a in axes.flat:
        a.axvline(release, color='#999999', ls='--', lw=1)
        a.set(xlabel='Simulation time (µs)', xlim=(0, t[-1]))
        a.grid(alpha=.16)
    fig.suptitle('HOPS single-pulse response · 1 ns native archives\nPulse release at 0.2 µs (dashed line) · collisionless development model', fontsize=13)
    fig.savefig(output/'relaxation.png', dpi=170)
    fig.savefig(output/'relaxation.pdf')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    measured = measure(args.run, args.output, args.reference)
    print(json.dumps({key: measured[key] for key in ('states', 'at_release', 'final', 'response', 'voltage_fall_tau_ns', 'native_checks')}, indent=2))
