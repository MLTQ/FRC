"""Free expansion of moving neutral, charged and electron-ion species volumes."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV, epsilon_0 as EPS0
from scipy.integrate import quad
from plasma_bench.atomic_frame import FixedFrameAtomicMoments
from plasma_bench.atomic_network import AtomicNetwork
from plasma_bench.atomic_moving import MovingAtomicElectrostatic
from plasma_bench.radiative_recombination import RadiativeRecombination
from .atomic_moments import ConstantIonization
from .implicit_plasma_records import provenance, freeze


def charged_expansion_time(scale, pressure_acceleration, plasma2, gamma=5/3):
    """Independent energy-integral quadrature for s''=A s^-gamma+omega_pi²."""
    def speed(s):
        return np.sqrt(-2*pressure_acceleration/(gamma-1)*np.expm1((1-gamma)*np.log(s))+2*plasma2*(s-1))
    def integrand(y):
        return 2*y/speed(1+y*y) if y else np.sqrt(2/(pressure_acceleration+plasma2))
    return quad(integrand, 0., np.sqrt(scale-1), epsabs=1e-15, epsrel=1e-11)[0], speed(scale)


def reference_amounts(faces, source, kind, half_width, temperature, time_s, scale=1., scale_rate=0.):
    """Exact continuum cell integrals split at actual analytic heads/fronts."""
    gamma = 5/3; result = np.zeros((len(faces)-1, source.size)); nodes, weights = np.polynomial.legendre.leggauss(8)
    sound = np.sqrt(gamma*EV*temperature/source.mass[0])
    if kind == 'neutral':
        front = half_width+2*sound*time_s/(gamma-1)
        breaks = [-front, -half_width+sound*time_s, half_width-sound*time_s, front]
        species = 0
    else:
        front = half_width*scale; breaks = [-front, front]; species = 1
    for i, (lo, hi) in enumerate(zip(faces[:-1], faces[1:])):
        cuts = np.r_[lo, [x for x in breaks if lo < x < hi], hi]
        for left, right in zip(cuts[:-1], cuts[1:]):
            x = .5*(left+right)+.5*(right-left)*nodes; weight = .5*(right-left)*weights
            if kind == 'neutral':
                if time_s == 0:
                    number = (abs(x) < half_width).astype(float); pressure = temperature*number; velocity = np.zeros_like(x)
                else:
                    xi = (abs(x)-half_width)/time_s; core = xi <= -sound; fan = (xi > -sound)&(xi < 2*sound/(gamma-1))
                    ratio = np.zeros_like(x); ratio[core] = 1.
                    ratio[fan] = 2/(gamma+1)*(1-(gamma-1)*xi[fan]/(2*sound))
                    number = ratio**(2/(gamma-1)); pressure = temperature*ratio**(2*gamma/(gamma-1)); velocity = np.zeros_like(x)
                    velocity[fan] = np.sign(x[fan])*2/(gamma+1)*(sound+xi[fan])
            else:
                inside = abs(x) < front; number = inside.astype(float)/scale
                velocity = scale_rate*x/scale
                pressure = np.where(inside, temperature*scale**(-gamma)*(1-(x/front)**2), 0.)
            result[i, species] += weight@number
            result[i, source.momentum[species, 2]] += weight@(source.mass[species]/source.momentum_scale*number*velocity)
            result[i, source.energy[species]] += weight@(pressure/(gamma-1)+.5*source.mass[species]/EV*number*velocity**2)
    return result


def run(destination, mesh_motion='affine'):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh moving-species output required')
    rr = RadiativeRecombination(dict(nuclear_charge=1, transitions=[dict(initial_charge=1, kind='power_law',
        A_cm3_s=1e-12, reference_temperature_ev=2., exponent=0., minimum_temperature_ev=None, maximum_temperature_ev=None)]))
    network = AtomicNetwork(ConstantIonization(), rr, ionization_enabled=False, recombination_enabled=False, radiation=False)
    source = FixedFrameAtomicMoments(network, 39.948*atomic_mass, [0., 0., 0.])
    identity = provenance(['examples.atomic_moving']); out.mkdir(parents=True); freeze(out, identity)
    report = dict(terminal=False, source_sha256=identity['source_sha256'], cases=[], mesh_motion=mesh_motion,
        scope='Moving full-species Euler-Poisson volumes with exact represented exterior vacuum. Neutral rarefaction and warm charged affine-expansion references; nonlinear electron-ion expansion with no formation or kinetic-validity claim.')
    began = time.monotonic(); finals = {}

    def save():
        report['elapsed_seconds'] = time.monotonic()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')

    save()
    for kind in ('neutral', 'charged', 'plasma'):
        half_width = .01 if kind == 'plasma' else .025
        density = 1e14 if kind == 'plasma' else (1e12 if kind == 'charged' else 1e19)
        temperature = 1. if kind == 'charged' else .1
        scale = 1.5; scale_rate = 0.
        if kind == 'charged':
            plasma2 = EV**2*density/(EPS0*source.mass[1]); pressure_acceleration = 2*EV*temperature/(source.mass[1]*half_width**2)
            duration, scale_rate = charged_expansion_time(scale, pressure_acceleration, plasma2)
        elif kind == 'neutral':
            sound = np.sqrt((5/3)*EV*temperature/source.mass[0]); duration = half_width/(4*sound)
        else:
            duration = 20e-9
        for cells in (32, 64, 128):
            faces = np.linspace(-half_width, half_width, cells+1)
            model = MovingAtomicElectrostatic(source, faces, density, (-.2, .2),
                left_electric=-EV*density*half_width/EPS0 if kind == 'charged' else 0., maximum_plasma_phase=.04*32/cells,
                mesh_motion=mesh_motion)
            if kind in ('neutral', 'charged'):
                amounts = reference_amounts(faces, source, kind, half_width, temperature, 0.)
            else:
                initial = source.initial(network.initial([0., 1.], 2.), .2, np.zeros((3, 3)))
                amounts = np.diff(faces)[:, None]*initial
            t = 0.; times = [0.]; positions = [faces.copy()]; states = [amounts.copy()]; failure = None; max_work = 0.
            while t < duration:
                dt = min(model.timestep(faces, amounts), duration-t)
                try:
                    candidate = model.step(faces, amounts, t, dt, reactions=False)
                except (ValueError, RuntimeError) as error:
                    failure = dict(step=len(times), time_s=t, dt_s=dt, kind=type(error).__name__, message=str(error)); break
                faces, amounts = candidate; t = min(t+dt, duration)
                times.append(t); positions.append(faces.copy()); states.append(amounts.copy())
                max_work = max(max_work, abs(model.last_diagnostic['electric_work_residual_ev_m_s']))
            positions = np.asarray(positions); states = np.asarray(states)
            rows = [model.totals(x, q) for x, q in zip(positions, states)]
            budgets = {key: np.asarray([r[key] for r in rows]) for key in rows[0]}
            fields = np.asarray([model.electric(q) for q in states])
            case = dict(name=f'{kind}_n{cells}', kind=kind, cells=cells, reference_density_m3=density,
                initial_half_width_m=half_width, duration_s=duration, completed_time_s=t, steps=len(times)-1,
                success=failure is None, failure=failure, maximum_dt_s=float(np.diff(times).max(initial=0.)),
                final_fronts_m=faces[[0, -1]].tolist(), minimum_cell_width_m=float(np.diff(positions, axis=1).min()),
                maximum_cell_width_ratio=float(np.max(np.diff(positions, axis=1).max(axis=1)/np.diff(positions, axis=1).min(axis=1))),
                maximum_nucleus_error=float(abs(budgets['nuclei']-budgets['nuclei'][0]).max()),
                maximum_charge_error=float(abs(budgets['charge']-budgets['charge'][0]).max()),
                maximum_momentum_error_scaled=float(np.linalg.norm(budgets['momentum_scaled']-budgets['momentum_scaled'][0], axis=1).max()),
                maximum_energy_error_ev_m=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()),
                maximum_energy_error_over_initial_total=float(abs(budgets['total_energy_ev']-budgets['total_energy_ev'][0]).max()/budgets['total_energy_ev'][0]),
                maximum_stage_electric_work_residual_ev_m_s=max_work,
                limited_reconstruction_fraction=model.limited_reconstruction_count/model.reconstruction_count,
                final_peak_electric_v_m=float(abs(fields[-1]).max()),
                exterior_species_state_maximum=float(abs(model.sample(faces, amounts, np.array([-.2, faces[0]-1e-10, faces[-1]+1e-10, .2]))).max()))
            history = dict(time_s=np.asarray(times), faces_m=positions, amounts=states, electric_v_m=fields, **budgets)
            if kind in ('neutral', 'charged') and failure is None:
                exact = reference_amounts(faces, source, kind, half_width, temperature, t, scale, scale_rate)
                front = half_width*scale if kind == 'charged' else half_width+3*sound*t
                species = 1 if kind == 'charged' else 0
                outside = reference_amounts(np.array([-.2, faces[0]]), source, kind, half_width, temperature, t, scale, scale_rate)
                outside += reference_amounts(np.array([faces[-1], .2]), source, kind, half_width, temperature, t, scale, scale_rate)
                case.update(analytic_front_m=front, maximum_front_error_m=float(abs(faces[[0, -1]]-np.array([-front, front])).max()),
                    final_number_l1_over_initial_number=float((abs(amounts[:, species]-exact[:, species]).sum()+outside[0, species])/(2*half_width)),
                    final_momentum_l1_scaled=float(abs(amounts[:, source.momentum[species, 2]]-exact[:, source.momentum[species, 2]]).sum()),
                    final_energy_l1_ev_m=float(abs(amounts[:, source.energy[species]]-exact[:, source.energy[species]]).sum()))
                if kind == 'charged':
                    exact_field = EV*density/EPS0*np.clip(faces/scale, -half_width, half_width)
                    case['maximum_final_electric_error_over_boundary_field'] = float(abs(fields[-1]-exact_field).max()/abs(model.left_electric))
                    history['analytic_final_electric_v_m'] = exact_field
                history['analytic_final_amounts'] = exact
            np.savez_compressed(out/(case['name']+'.npz'), **history)
            report['cases'].append(case); finals[kind, cells] = (model, faces, amounts, history); save()
            print('Moving species', case, flush=True)
    grid = np.linspace(-.19, .19, 10001); comparisons = []
    for coarse, fine in ((32, 64), (64, 128)):
        a, b = (finals['plasma', k] for k in (coarse, fine))
        if any(not r['success'] for r in report['cases'] if r['kind'] == 'plasma' and r['cells'] in (coarse, fine)):
            continue
        da, db = (model.sample(faces, amounts, grid) for model, faces, amounts, _ in (a, b))
        ea, eb = (np.interp(grid, faces, model.electric(amounts)) for model, faces, amounts, _ in (a, b))
        comparisons.append(dict(coarse_cells=coarse, fine_cells=fine,
            final_number_l1_over_initial_number=np.trapezoid(abs(da[:, source.network.electrons]-db[:, source.network.electrons]), grid)/.02,
            final_electric_relative_l2=float(np.linalg.norm(ea-eb)/np.linalg.norm(eb)),
            final_front_difference_m=float(abs(a[1][[0, -1]]-b[1][[0, -1]]).max())))
    report['plasma_spatial_sensitivity'] = comparisons
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for cells in (32, 64, 128):
        for kind, ax in (('neutral', axes[0, 0]), ('charged', axes[0, 1])):
            model, x, q, h = finals[kind, cells]; species = 0 if kind == 'neutral' else 1
            ax.plot(.5*(x[:-1]+x[1:]), q[:, species]/np.diff(x), label=f'{cells} cells')
            if cells == 128 and 'analytic_final_amounts' in h:
                ax.plot(.5*(x[:-1]+x[1:]), h['analytic_final_amounts'][:, species]/np.diff(x), 'k--', label='Independent continuum reference')
        model, x, q, h = finals['plasma', cells]
        completed = next(r['success'] for r in report['cases'] if r['kind'] == 'plasma' and r['cells'] == cells)
        status = 'complete' if completed else 'stopped'
        axes[1, 0].plot(x, model.electric(q), label=f'{cells} cells: {h["time_s"][-1]*1e9:.3f} ns ({status})')
        axes[1, 1].plot(h['time_s']*1e9, (h['total_energy_ev']-h['total_energy_ev'][0])/h['total_energy_ev'][0], label=f'{cells} cells ({status})')
    labels = [('Neutral expansion', 'Number / Nref'), ('Charged affine expansion', 'Number / Nref'),
              ('Last accepted plasma fields: different times', 'Electric field [V/m]'), ('Electron–ion energy budget', 'Relative total-energy drift')]
    for i, (ax, (title, ylabel)) in enumerate(zip(axes.flat, labels)):
        ax.set(title=title, xlabel='Time [ns]' if i == 3 else 'Axial coordinate [m]', ylabel=ylabel); ax.legend(fontsize=7); ax.grid(alpha=.2)
    fig.suptitle('Moving species and electric field with exact exterior vacuum\nPlanar fluid model; front error and conservation measured separately')
    fig.savefig(out/'moving_species.png', dpi=170); fig.savefig(out/'moving_species.pdf'); plt.close(fig)
    if provenance(['examples.atomic_moving']) != identity:
        raise RuntimeError('Source changed during moving species evolution')
    report['terminal'] = True; report['success'] = all(r['success'] for r in report['cases']); save(); return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('destination'); p.add_argument('--mesh-motion', choices=('affine', 'barycentric'), default='affine'); a = p.parse_args()
    r = run(a.destination, a.mesh_motion); raise SystemExit(0 if r['success'] else 1)
