"""Run one configured biased-coil experiment with the complete material solver."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME
from plasma_bench.material_experiment import build
from plasma_bench.material_carrier_current import MaterialCurrentState
from plasma_bench.coulomb_transport import regime
from .implicit_plasma_records import provenance, freeze, sha


def run(config_path, destination):
    config = json.loads(Path(config_path).read_text()); out = Path(destination)
    if out.exists():
        raise ValueError('A fresh experiment output directory is required')
    integrator = config.get('integrator', 'rk4')
    identity = provenance(['examples.material_experiment']); started = time.monotonic()
    out.mkdir(parents=True); freeze(out, identity)
    (out/'config.json').write_text(json.dumps(config, indent=2)+'\n')
    record = dict(terminal=False, success=False, completed_steps=0, config=config,
        source_sha256=identity['source_sha256'], integrator=integrator, records=[],
        scope='Complete moving-plasma equations with an assumed pre-ionized inventory and prescribed biased coils. Preparation, ionization, neutrals, physical edge/transport validity and convergence are not established. No formation or laboratory recommendation follows from this exploratory run.')

    def save():
        record.update(elapsed_seconds=time.monotonic()-started,
                      peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        (out/'summary.json').write_text(json.dumps(record, indent=2)+'\n')

    save()
    try:
        fluid, model, metadata = build(config); record['resolved_inputs'] = metadata
        steps = config['evolution']['steps']; duration = config['evolution']['duration_s']; dt = duration/steps
        state = MaterialCurrentState(fluid.initial, np.zeros(model.space.magnetic_count))
        previous = model.response(state, time=0.); energy0 = previous['total_energy']
        work = 0.; response_work = 0.; heat = 0.
        group_work = np.zeros(len(model.source_groups)); group_response_work = group_work.copy()
        field_per_amp = np.array([group['vacuum_center_field_per_amp_t'] for group in metadata['source_groups']])
        record.update(steps=steps, dt=dt, duration=duration, initial_total_energy_j=energy0)

        def observe(s, t, a):
            f = a['fields']; n = f['n'].ravel(); temperature = f['temperature'].reshape(-1, 2)
            j = E*n[:, None]*a['drift']; b = a['magnetic_field']; vol = a['volume']
            points = f['points'].reshape(-1, 2)
            source_b = a['applied_magnetic_field']
            psi = model.space.flux@s.magnetic+points[:, 0]*a['applied_vector_potential'][:, 1]
            drift_ratio = np.linalg.norm(a['drift'], axis=1)/np.sqrt(E*temperature[:, 0]/ME)
            k = int(drift_ratio.argmax()); diagnostics = config['diagnostic_scales']
            local = regime(n[k], *temperature[k], fluid.ion_mass,
                field=np.linalg.norm(b[k]), current_density=np.linalg.norm(j[k]),
                length=diagnostics['reference_length_m'], timescale=diagnostics['reference_time_s'])
            boundary = s.fluid.coordinates[fluid.domain.boundary]
            return dict(time=t, source_current_a=a['source_current_a'],
                nominal_vacuum_center_field_t=float(a['source_group_currents_a']@field_per_amp),
                source_group_currents_a=a['source_group_currents_a'].tolist(),
                source_group_center_fields_t=(a['source_group_currents_a']*field_per_amp).tolist(),
                source_group_linkages_wb=a['source_group_linkages_wb'].tolist(),
                source_group_external_work_j=group_work.tolist(), source_group_response_work_j=group_response_work.tolist() if integrator=='rk4' else None,
                source_linkage_wb=a['source_linkage_wb'], total_energy_j=a['total_energy'],
                response_energy_j=a['response_energy'], external_work_j=work,
                source_response_work_j=response_work if integrator=='rk4' else None, joule_heat_j=heat,
                relative_external_work_residual=(a['total_energy']-energy0-work)/energy0,
                relative_response_work_residual=(a['response_energy']-energy0-response_work)/energy0 if integrator=='rk4' else None,
                induced_magnetic_energy_j=a['magnetic_energy'], current_energy_j=a['current_energy'],
                peak_current_components_a_m2=abs(j).max(axis=0).tolist(),
                peak_field_components_t=abs(b).max(axis=0).tolist(),
                plasma_bz_range_t=[float(b[:, 2].min()), float(b[:, 2].max())],
                plasma_flux_range_wb_per_radian=[float(psi.min()), float(psi.max())],
                plasma_volume_fraction_bz_opposing_local_source=float(np.sum(vol*(b[:, 2]*source_b[:, 2]<0))/vol.sum()),
                particle_mean_temperature_ev=(np.sum(fluid.particles[..., None]*f['temperature'], axis=(0, 1))/fluid.particles.sum()).tolist(),
                minimum_temperature_ev=temperature.min(axis=0).tolist(),
                maximum_drift_over_electron_thermal_speed=float(drift_ratio[k]),
                local_regime_at_max_drift_ratio=dict(position_r_z_m=points[k].tolist(), **local),
                maximum_plasma_radius_m=float(boundary[:, 0].max()),
                midplane_plasma_radius_m=float(np.interp(0., model.space.mesh.zf, s.fluid.coordinates[fluid.domain.ids[-1], 0])),
                particle_rms_radius_m=float(np.sqrt(np.sum(fluid.particles*f['points'][..., 0]**2)/fluid.particles.sum())),
                maximum_plasma_half_length_m=float(abs(boundary[:, 1]).max()),
                minimum_specific_entropy=s.fluid.entropy.min(axis=(0, 1)).tolist(),
                minimum_vacuum_map_jacobian=float(a['vacuum_geometry']['determinant'].min()),
                maximum_transport_outflow_s=a['maximum_transport_outflow_s'],
                exchange_coulomb_log_range=a['exchange_coulomb_log_range'],
                pair_power_residual_w=a['pair_power_residual'], residuals=a['residuals'],
                split_diagnostics=getattr(model, 'split_diagnostics', None))

        def arrays(s):
            return dict(magnetic=s.magnetic, **{k: getattr(s.fluid, k) for k in s.fluid.__dataclass_fields__})

        np.savez_compressed(out/'initial.npz', **arrays(state), particles=fluid.particles, ids=fluid.ids)
        record['initial'] = observe(state, 0., previous); save(); stride = max(1, steps//80)
        for step in range(1, steps+1):
            next_state = model.step(state, dt, (step-1)*dt)
            a = model.response(next_state, time=step*dt)
            if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 4*1024**3 or shutil.disk_usage('.').free < 256*1024**2:
                raise RuntimeError('Own4GiB/256MiB resource guard')
            work += .5*dt*(previous['external_power']+a['external_power'])
            if integrator == 'rk4':
                response_work += .5*dt*(previous['source_response_power']+a['source_response_power'])
            group_work += .5*dt*(previous['source_group_external_powers_w']+a['source_group_external_powers_w'])
            if integrator == 'rk4':
                group_response_work += .5*dt*(previous['source_group_response_powers_w']+a['source_group_response_powers_w'])
            heat += .5*dt*(previous['joule_power']+a['joule_power'])
            state = next_state; previous = a; record['completed_steps'] = step
            raw = out/f'step_{step:05}.npz'; np.savez_compressed(raw, **arrays(state))
            if step % stride == 0 or step in (1, steps):
                row = observe(state, step*dt, a)
                row.update(step=step, arrays_file=raw.name, arrays_sha256=sha(raw))
                record['records'].append(row); save()
            if step == 1 or step % max(1, steps//8) == 0:
                print('material experiment', step, '/', steps, flush=True)
        record['success'] = True
    except Exception as error:
        record.update(failure_type=type(error).__name__, failure_message=str(error), failed_step=record['completed_steps']+1)
    record['terminal'] = True
    if provenance(['examples.material_experiment']) != identity:
        raise RuntimeError('Source changed during experiment')
    save(); print('material experiment terminal', record['success'], record.get('failure_message', ''), flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('config'); parser.add_argument('destination')
    args = parser.parse_args(); result = run(args.config, args.destination)
    raise SystemExit(0 if result['success'] else 1)
