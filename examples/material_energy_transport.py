"""Complete graph heat evolutions against independent continuum solutions."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from plasma_bench.material_energy_diffusion import diffuse_energy
from .implicit_plasma_records import provenance, freeze


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('A fresh output directory is required')
    identity = provenance(['examples.material_energy_transport']); started = time.monotonic()
    out.mkdir(parents=True); freeze(out, identity)
    record = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[],
        scope='Independent heat-equation and two-reservoir limits of the numerical thermal mixing operator. Not a physical plasma heat closure or full coupled-solver validation.')
    times = np.linspace(0., .5, 51); gamma = 5/3
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for mode in ('fixed_diffusivity', 'fixed_drift'):
        for count in (16, 32, 64, 128):
            dx = 1/count; x = (np.arange(count)+.5)*dx
            diffusivity = .1 if mode == 'fixed_diffusivity' else .5*dx
            mass = np.full(count, dx); tail = np.arange(count); head = (tail+1)%count
            conductance = np.full(count, diffusivity/dx)
            initial = 1+.8*np.sinc(dx)*np.cos(2*np.pi*x)
            history = [initial]
            for dt in np.diff(times):
                history.append(diffuse_energy(mass, tail, head, conductance, history[-1], dt))
            history = np.array(history)
            exact = 1+.8*np.sinc(dx)*np.exp(-4*np.pi**2*diffusivity*times[:, None])*np.cos(2*np.pi*x)
            energy = history@mass; entropy = np.log(history)@mass/(gamma-1)
            error = np.sqrt(np.sum(mass*(history-exact)**2, axis=1))
            initial_distance = np.sqrt(np.sum(mass*(history[-1]-initial)**2))
            row = dict(mode=mode, cells=count, diffusivity=diffusivity,
                final_l2_error=float(error[-1]), maximum_l2_error=float(error.max()),
                final_distance_from_initial=float(initial_distance),
                maximum_relative_energy_change=float(np.max(abs(energy/energy[0]-1))),
                minimum_energy_per_particle=float(history.min()),
                total_entropy_change=float(entropy[-1]-entropy[0]),
                minimum_sampled_entropy_increment=float(np.diff(entropy).min()))
            record['cases'].append(row)
            np.savez_compressed(out/f'{mode}_{count}.npz', x=x, times=times,
                internal=history, exact=exact, particle_counts=mass, entropy=entropy, total_energy=energy)
            column = 0 if mode == 'fixed_diffusivity' else 1
            axes[0, column].plot(x, history[-1], label=f'{count} cells')
            if count == 128:
                axes[0, column].plot(x, exact[-1], 'k--', label='Continuum heat equation')
                if mode == 'fixed_drift':
                    axes[0, column].plot(x, initial, color='.5', ls=':', label='Zero-diffusion limit')
        rows = [r for r in record['cases'] if r['mode'] == mode]
        axes[1, column].loglog([r['cells'] for r in rows], [r['final_l2_error'] for r in rows], 'o-', label='Continuum discretization error')
        if mode == 'fixed_drift':
            axes[1, column].loglog([r['cells'] for r in rows], [r['final_distance_from_initial'] for r in rows], 's-', label='Change from initial profile')
    mass = np.array([1., 3.]); conductance = np.array([.7]); h = np.array([4., .2])
    mean = float(mass@h/mass.sum()); rate = .7*np.sum(1/mass)
    history = [h]
    for dt in np.diff(times):
        history.append(diffuse_energy(mass, [0], [1], conductance, history[-1], dt))
    history = np.array(history); exact = mean+(h-mean)*np.exp(-rate*times[:, None])
    record['unequal_reservoirs'] = dict(particles=mass.tolist(), initial_internal=h.tolist(),
        conductance=.7, analytic_decay_rate=rate,
        maximum_absolute_error=float(abs(history-exact).max()),
        maximum_relative_energy_change=float(np.max(abs((history@mass)/(h@mass)-1))),
        minimum_sampled_entropy_increment=float(np.diff(np.log(history)@mass/(gamma-1)).min()))
    np.savez_compressed(out/'unequal_reservoirs.npz', times=times, internal=history, exact=exact, particle_counts=mass)
    for column, title in enumerate(('Fixed D = 0.1', 'Upwind-scale D = |u| dx / 2, |u| = 1')):
        axes[0, column].set(title=title, xlabel='Periodic position', ylabel='Internal energy per particle')
        axes[1, column].set(xlabel='Cells', ylabel='Particle-weighted L2 difference')
    for ax in axes.flat:
        ax.legend(frameon=False); ax.grid(alpha=.2)
    fig.suptitle('Conservative numerical thermal mixing: independent diffusion limits\nNo physical plasma conductivity is inferred')
    fig.savefig(out/'thermal_mixing_limits.png', dpi=170); plt.close(fig)
    record.update(terminal=True, elapsed_seconds=time.monotonic()-started)
    (out/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    print('Thermal transport limits complete', out, flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    run(parser.parse_args().destination)
