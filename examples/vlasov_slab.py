"""Full finite-current kinetic slab evolution and independent ballistic endpoint."""
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from scipy.constants import atomic_mass, electron_mass, elementary_charge as EV
from scipy.special import erfc
from plasma_bench.vlasov_slab import KineticSpecies, VlasovSlab
from plasma_bench.ballistic_slab import cell_moments
from .implicit_plasma_records import provenance, freeze


def snapshot(path, model, state):
    data = dict(x_edges_m=model.edges, electric_edges_v_m=model.electric_field(state))
    for s, f in zip(model.species, state):
        data[s.name+'_g'] = f; data[s.name+'_s_edges'] = s.s_edges
        data.update({s.name+'_'+k: v for k, v in s.profiles(f, model.density).items()})
    np.savez_compressed(path, **data)


def run(destination, cells, velocity_cells, self_consistent, sigma=8., cfl=.65):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh output directory required')
    if shutil.disk_usage(out.parent).free < 384*1024**2:
        raise OSError('Insufficient storage for output and reserved working space')
    identity = provenance(['examples.vlasov_slab'])
    out.mkdir(); freeze(out, identity)
    species = [KineticSpecies('ion', 39.948*atomic_mass-electron_mass, EV, .2, velocity_cells, sigma),
               KineticSpecies('electron', electron_mass, -EV, 2., velocity_cells, sigma)]
    model = VlasovSlab(species, cells, self_consistent=self_consistent)
    state = model.initial(); snapshot(out/'initial.npz', model, state)
    total, energy, field_energy = model.budget(state)
    initial_energy, initial_numbers = energy, total[:, 0].copy()
    target = 20e-9
    tail_probability = erfc(sigma/np.sqrt(2))
    tail_second = tail_probability+2*sigma*np.exp(-sigma*sigma/2)/np.sqrt(2*np.pi)
    exact_initial_energy = 2*model.half_width*model.density*sum(
        EV*s.temperature_ev*(1.5-tail_probability-.5*tail_second) for s in species)
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
        self_consistent=self_consistent, spatial_cells=cells, velocity_cells=velocity_cells,
        cfl=cfl, velocity_sigma=sigma, gaussian_omitted_probability=float(tail_probability),
        gaussian_omitted_second_moment=float(tail_second), requested_time_s=target,
        initial_p0_energy_relative_to_exact_truncated=float(energy/exact_initial_energy-1),
        density_m3=model.density, half_width_m=model.half_width, ends_m=[-.2, .2],
        species=[dict(name=s.name, mass_kg=s.mass, charge_c=s.charge, temperature_ev=s.temperature_ev,
                      thermal_speed_m_s=s.thermal_speed) for s in species],
        scope='Collisionless electrostatic slab; free case is a mathematical limit; no formation or lab claim.')
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    now, dt_history, history, profile_times, profiles = 0., [], [], [], []
    loss_x, loss_v, impulse, failure = np.zeros((2, 3)), np.zeros((2, 3)), 0., None
    started = time.perf_counter()
    def record():
        totals, e, ef = model.budget(state)
        history.append(np.r_[now, e, ef, totals.ravel(), loss_x.ravel(), loss_v.ravel(), impulse])
    def profile():
        profile_times.append(now)
        item = {'electric_edges_v_m': model.electric_field(state)}
        for s, f in zip(species, state):
            item.update({s.name+'_'+k: v for k, v in s.profiles(f, model.density).items()})
        profiles.append(item)
    record(); profile()
    while now < target:
        try:
            next_profile = min(len(profile_times)*1e-9, target)
            dt = min(model.timestep(state, cfl=cfl), next_profile-now, target-now)
            if dt <= 0 or not np.isfinite(dt) or now+dt == now:
                raise ValueError('Kinetic clock cannot advance')
            candidate, xlost, vlost, boundary_impulse = model.step(state, dt)
        except (ValueError, ArithmeticError, RuntimeError) as exc:
            failure = dict(time_s=now, message=str(exc)); break
        state = candidate; now = min(now+dt, next_profile, target)
        loss_x += xlost; loss_v += vlost; impulse += boundary_impulse
        dt_history.append(dt); record()
        if now == next_profile:
            profile()
        if len(dt_history) % 250 == 0:
            print('kinetic', self_consistent, cells, velocity_cells, len(dt_history), now, flush=True)
            if shutil.disk_usage(out).free < 320*1024**2:
                failure = dict(time_s=now, message='Storage reserve reached; accepted prefix saved'); break
    snapshot(out/'final.npz', model, state)
    rows = np.asarray(history)
    totals, xlosses, vlosses = rows[:, 3:9].reshape(-1, 2, 3), rows[:, 9:15].reshape(-1, 2, 3), rows[:, 15:21].reshape(-1, 2, 3)
    balanced = totals+xlosses+vlosses
    corrected_energy = rows[:, 1]+(xlosses+vlosses)[:, :, 2].sum(axis=1)
    corrected_momentum = balanced[:, :, 1].sum(axis=1)-rows[:, 21]
    np.savez_compressed(out/'history.npz', time_s=rows[:, 0], dt_s=dt_history,
        energy_j_m2=rows[:, 1], electric_energy_j_m2=rows[:, 2], species_npe=totals,
        spatial_outflow_npe=xlosses, velocity_outflow_npe=vlosses,
        electric_boundary_impulse_kg_m_s=rows[:, 21])
    np.savez_compressed(out/'profiles.npz', time_s=profile_times, x_edges_m=model.edges,
        **{k: np.stack([p[k] for p in profiles]) for k in profiles[0]})
    report.update(terminal=True, success=bool(failure is None and now == target), failure=failure,
        completed_time_s=now, steps=len(dt_history), elapsed_seconds=time.perf_counter()-started,
        maximum_relative_energy_balance_error=float(np.max(abs(corrected_energy-initial_energy))/initial_energy),
        maximum_relative_number_balance_error=float(np.max(abs(balanced[:, :, 0]-initial_numbers)/initial_numbers)),
        maximum_momentum_balance_error_kg_m_s=float(np.max(abs(corrected_momentum-corrected_momentum[0]))),
        final_spatial_outflow_npe=loss_x.tolist(), final_velocity_outflow_npe=loss_v.tolist(),
        maximum_field_v_m=float(np.max(abs(model.electric_field(state)))))
    if report['success'] and not self_consistent:
        exact, metrics = {}, {}
        for s, f in zip(species, state):
            reference = cell_moments(model.edges, model.half_width, s.thermal_speed, now, sigma)
            actual = s.raw(f)
            exact[s.name+'_raw'] = reference
            norm = np.sqrt(np.sum(reference**2, axis=0)*model.dx)
            errors = np.sqrt(np.sum((actual-reference)**2, axis=0)*model.dx)/norm
            metrics[s.name] = dict(raw_moment_relative_l2=errors.tolist(),
                reference_total_number_relative_error=float(reference[:, 0].sum()*model.dx/(2*model.half_width*(1-tail_probability))-1))
        np.savez_compressed(out/'exact_characteristics.npz', x_edges_m=model.edges, **exact)
        report['independent_ballistic'] = metrics
    if identity != provenance(['examples.vlasov_slab']):
        report['source_changed_during_run'] = True
        report['success'] = False
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(report, flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination')
    parser.add_argument('--cells', type=int, required=True)
    parser.add_argument('--velocity-cells', type=int, required=True)
    parser.add_argument('--free', action='store_true'); parser.add_argument('--sigma', type=float, default=8.)
    parser.add_argument('--cfl', type=float, default=.65)
    args = parser.parse_args()
    result = run(args.destination, args.cells, args.velocity_cells, not args.free, args.sigma, args.cfl)
    raise SystemExit(0 if result['success'] else 1)
