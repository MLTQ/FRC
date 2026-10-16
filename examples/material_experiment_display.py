"""Show completed configured coil experiments, including the actual exterior B."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.constants import elementary_charge as E
from plasma_bench.material_experiment import build as experiment, source_configs
from plasma_bench.material_plasma_samples import material_samples
from .material_coil_plot import state_at
from .material_vector_display import outline, trace
from .implicit_plasma_records import sha


def circle(radius, z):
    angle = np.linspace(0, 2*np.pi, 97)
    return np.column_stack((np.full_like(angle, z), radius*np.cos(angle), radius*np.sin(angle)))


def build(source, destination, *, allow_partial=False, frame_interval_ns=None):
    source = Path(source); out = Path(destination); d = json.loads((source/'summary.json').read_text())
    if out.exists() or not d['terminal'] or (not d['success'] and not allow_partial) or d['completed_steps'] < 1:
        raise ValueError('Terminal experiment and fresh display required; failed prefixes need --allow-partial')
    last_step = d['completed_steps']
    status = 'prescribed biased coils · exploratory plasma model' if d['success'] else f"Stopped at step {last_step+1}: {d.get('failure_message', 'simulation failed')}. Accepted states only."
    config = d['config']; fluid, model, metadata = experiment(config)
    tube = config['tube']; plasma = config['plasma']; grid = config['mesh']
    viscosity = config.get('numerical_bulk_viscosity')
    viscosity_note = (' No numerical bulk viscosity.' if viscosity is None else
        f" Numerical compression pressure is enabled (linear {viscosity['linear']:g}, quadratic {viscosity['quadratic']:g}); its ion heating is numerical, not a physical transport model.")
    external_frames = frame_interval_ns is not None
    if external_frames:
        interval = frame_interval_ns*1e-9
        if not np.isfinite(interval) or interval < d['dt']:
            raise ValueError('Display interval must be finite and at least the saved integration step')
        targets = np.arange(int(np.floor(last_step*d['dt']/interval+1e-10))+1)*interval
        steps = set(np.clip(np.rint(targets/d['dt']), 0, last_step).astype(int)) | {last_step}
    else:
        steps = set(np.round(np.linspace(0, last_step, 8)).astype(int))
        for source_config in source_configs(config):
            pulse = source_config['pulse_center_t']
            if pulse['kind'] == 'Unipolar':
                for t in (pulse['start'], pulse['start']+pulse['rise'], pulse['start']+pulse['rise']+pulse['hold'], pulse['start']+pulse['rise']+pulse['hold']+pulse['fall']):
                    steps.add(min(last_step, max(0, round(t/d['dt']))))
    out.mkdir(parents=True)
    if external_frames: (out/'frames').mkdir()
    extent = max(tube['length_m']/2, max(c[0] for c in model.coils)/.6)*1.05
    def normalized(frame):
        cloud = np.asarray(frame['current']).copy(); cloud[:, :3] /= extent
        return dict(frame, lines=[dict(line, points=(np.asarray(line['points'])/extent).tolist()) for line in frame['lines']],
                    outline=[(np.asarray(p)/extent).tolist() for p in frame['outline']], current=cloud.tolist())
    selected = None
    frames = []; inputs = {str(source/'summary.json'): sha(source/'summary.json')}
    q = grid['quadrature']; local = np.array([q+1, 2*q-2, (q-2)*q+1, (q-1)*q-2])
    indices = (np.arange(grid['nr']*grid['nz'])[:, None]*q*q+local[None, :]).ravel()
    coils = [circle(radius, z) for radius, z, turns in model.coils]
    coil_groups = [dict(label=g['name'], color=('#db9160', '#96b5ef', '#c6a5d8')[i % 3],
                       paths=[circle(radius, z) for radius, z, turns in g['coils_radius_z_turns']])
                   for i, g in enumerate(metadata['source_groups'])]
    chamber = [circle(tube['radius_m'], z) for z in (-tube['length_m']/2, 0., tube['length_m']/2)]
    for angle in np.arange(8)*np.pi/4:
        chamber.append(np.array([[-tube['length_m']/2, tube['radius_m']*np.cos(angle), tube['radius_m']*np.sin(angle)],
                                  [tube['length_m']/2, tube['radius_m']*np.cos(angle), tube['radius_m']*np.sin(angle)]]))
    for step in sorted(steps):
        raw = source/('initial.npz' if step == 0 else f'step_{step:05}.npz'); inputs[str(raw)] = sha(raw)
        state = state_at(source, step); t = step*d['dt']; a = model.response(state, rate=False, time=t)
        points = a['fields']['points'].reshape(-1, 2); current = E*a['fields']['n'].ravel()[:, None]*a['drift']
        vacuum = a['vacuum_response']; vp = vacuum['points']
        keep = (vp[:, 0] < 1.5*tube['radius_m']) & (abs(vp[:, 1]) < tube['length_m']/2)
        vp = vp[keep]; vb = vacuum['magnetic_field'][keep]
        source_field = model.applied_source(state.fluid.coordinates, vp, t)['magnetic']
        vb = np.column_stack((vb[:, 0], np.zeros(len(vb)), vb[:, 1]))+source_field
        all_points = np.concatenate((points, vp)); all_b = np.concatenate((a['magnetic_field'], vb))
        paths = []
        for sign in (1., -1.):
            paths.extend(trace(all_points, sign*all_b, 1.6*plasma['radius_m'], plasma['half_length_m']))
        colors = np.column_stack((np.linalg.norm(current, axis=1), current))
        selected_samples = np.unique(np.r_[indices, np.argmax(abs(colors), axis=0)]); cloud = []
        for angle in np.arange(12)*np.pi/6:
            cloud.extend(np.column_stack((points[selected_samples, 1], points[selected_samples, 0]*np.cos(angle),
                         points[selected_samples, 0]*np.sin(angle), colors[selected_samples])).tolist())
        peak_line_b = max((max(p['field']) for p in paths), default=0.)
        frame = dict(step=int(step), time_us=t*1e6, lines=paths, current=cloud,
            outline=[p.tolist() for p in outline(state.fluid.coordinates[fluid.domain.ids])],
            bmax=max(peak_line_b, float(np.linalg.norm(a['magnetic_field'], axis=1).max())),
            jmax=abs(colors).max(axis=0).tolist(),
            drive_center_t=float(sum(current*g['vacuum_center_field_per_amp_t'] for current, g in zip(a['source_group_currents_a'], metadata['source_groups']))),
            drive_groups=[dict(name=g['name'], center_t=float(current*g['vacuum_center_field_per_amp_t'])) for current, g in zip(a['source_group_currents_a'], metadata['source_groups'])],
            midplane_radius_m=float(np.interp(0., model.space.mesh.zf, state.fluid.coordinates[fluid.domain.ids[-1], 0])),
            maximum_radius_m=float(state.fluid.coordinates[fluid.domain.boundary, 0].max()),
            temperature=(np.sum(fluid.particles[..., None]*a['fields']['temperature'], axis=(0, 1))/fluid.particles.sum()).tolist())
        if selected is None or frame['jmax'][0] > selected['jmax'][0]: selected = frame
        if external_frames:
            filename = f'frames/step_{step:05}.json'
            (out/filename).write_text(json.dumps(normalized(frame), separators=(',', ':'))+'\n')
            frames.append(dict({k: v for k, v in frame.items() if k not in ('lines', 'current', 'outline')}, url=filename))
        else:
            frames.append(frame)
        print('experiment display state', step, flush=True)
    # Choose the largest-current saved state for the scientific still. The HTML
    # still exposes the initial bias, entire pulse and final relaxation separately.
    fig = plt.figure(figsize=(15, 7), constrained_layout=True)
    axes = [fig.add_subplot(1, 2, k+1, projection='3d') for k in range(2)]
    norm = Normalize(0., selected['bmax']*1e3)
    for path in selected['lines']:
        p = np.asarray(path['points']); segments = np.stack((p[:-1], p[1:]), axis=1)
        artist = Line3DCollection(segments, cmap='viridis', norm=norm, linewidths=.7)
        artist.set_array(np.asarray(path['field'][1:])*1e3); axes[0].add_collection3d(artist)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'), ax=axes[0], shrink=.6, label='|B| [mT]')
    cloud = np.asarray(selected['current'])
    artist = axes[1].scatter(*cloud[:, :3].T, c=cloud[:, 3]/1000, s=3, alpha=.6,
                            cmap='viridis', vmin=0., vmax=selected['jmax'][0]/1000)
    fig.colorbar(artist, ax=axes[1], shrink=.6, label='|J| [kA/m²]')
    for ax, title in zip(axes, ('Plasma and vacuum magnetic field lines', 'Actual plasma current density')):
        for paths, color, width in ((selected['outline'], '#5d7c83', .3), (chamber, '#738088', .4)):
            for path in paths: ax.plot(*np.asarray(path).T, color=color, lw=width, alpha=.7)
        for group in coil_groups:
            for i, path in enumerate(group['paths']):
                ax.plot(*path.T, color=group['color'], lw=1.2, alpha=.8, label=group['label'] if i == 0 else None)
        ax.legend(loc='upper left', fontsize=8, frameon=False)
        ax.set(xlim=(-extent, extent), ylim=(-.3*extent, .3*extent), zlim=(-.3*extent, .3*extent),
               xlabel='Axial z [m]', ylabel='x [m]', zlabel='y [m]', title=title)
        ax.set_box_aspect((2, .6, .6)); ax.view_init(22, -67)
    fig.suptitle(f"Biased coil pulse · actual saved state at {selected['time_us']:.3f} µs\n{tube['length_m']:g} m tube · assumed pre-ionized {plasma['gas']} · coarse exploratory calculation", fontsize=15)
    fig.savefig(out/'material_experiment_3d.png', dpi=180); fig.savefig(out/'material_experiment_3d.pdf'); plt.close(fig)
    # Show signed axial field and azimuthal current directly in the moving
    # meridional plane. Contours are actual sampled flux, not a fitted FRC shape.
    fig, axes = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    show = [frames[len(frames)//3], frames[2*len(frames)//3], frames[-1]]
    rr, zz = np.meshgrid(np.linspace(plasma['radius_m']*1e-5, plasma['radius_m']*(1-1e-5), 65),
                         np.linspace(-plasma['half_length_m']*(1-1e-5), plasma['half_length_m']*(1-1e-5), 257), indexing='ij')
    n0 = metadata['central_density_m3']*((1-(rr/plasma['radius_m'])**2)*(1-(zz/plasma['half_length_m'])**2))**1.5
    for column, frame in enumerate(show):
        state = state_at(source, frame['step'])
        samples = material_samples(model, state, rr, zz, n0, time=frame['step']*d['dt'])
        positions = samples['points']; psi = samples['flux']
        for row, (values, label) in enumerate(((samples['magnetic'][..., 2]*1e3, 'Bz [mT]'), (samples['current'][..., 1]/1000, 'Jθ [kA/m²]'))):
            ax = axes[row, column]; cap = float(abs(values).max())
            artist = ax.pcolormesh(positions[..., 1], positions[..., 0], values, shading='nearest', cmap='RdBu_r', vmin=-cap, vmax=cap, rasterized=True)
            if row == 0 and np.ptp(psi) > 0:
                ax.contour(positions[..., 1], positions[..., 0], psi, levels=np.linspace(psi.min(), psi.max(), 17)[1:-1], colors='#243a42', linewidths=.55)
                if psi.min() < 0 < psi.max(): ax.contour(positions[..., 1], positions[..., 0], psi, levels=[0.], colors='#f2c65c', linewidths=1.4)
            ax.set(xlabel='Axial z [m]', ylabel='Radius r [m]', title=f"{frame['time_us']:.3f} µs · {label}", aspect='auto')
            fig.colorbar(artist, ax=ax, shrink=.6, label=label)
    fig.suptitle(f"Actual biased-pulse fields on the moving plasma\n{grid['nr']}×{grid['nz']} simulation; lines: poloidal flux; radial scale expanded for readability", fontsize=15)
    fig.savefig(out/'material_experiment_sections.png', dpi=180); fig.savefig(out/'material_experiment_sections.pdf'); plt.close(fig)
    if not external_frames: frames = [normalized(frame) for frame in frames]
    data = dict(solver_dt_s=d['dt'], solver_steps=d['steps'], solver_completed_steps=d['completed_steps'], frames=frames, bmax=max(f['bmax'] for f in frames), jmax=np.max([f['jmax'] for f in frames], axis=0).tolist(),
        display_interval_ns=frame_interval_ns, initial_frame=0 if external_frames else len(frames)-1,
        playback_fps=10 if external_frames else 1,
        chamber=[(p/extent).tolist() for p in chamber], coils=[(p/extent).tolist() for p in coils], physical_scale_m=extent, initial_zoom=1.5,
        coil_groups=[dict(label=g['label'], color=g['color'], paths=[(p/extent).tolist() for p in g['paths']]) for g in coil_groups],
        status_note=status,
        scope=f"Source: {source.name}. {last_step} accepted of {d['steps']} requested coupled steps. {status}. Gas {plasma['gas']}, neutral preparation inventory {plasma['neutral_fill_mtorr']} mTorr at {plasma['neutral_temperature_k']} K, all assumed prepared as a compact singly ionized plasma. Central density {metadata['central_density_m3']:.4g}/m³; central Te/Ti {plasma['central_temperature_ev']} eV. Grid {grid['nr']}×{grid['nz']}, degree {grid['degree']}. Actual moving plasma and magnetic-only vacuum; prescribed resistance, conditional Coulomb exchange, explicit entropy diffusion.{viscosity_note} Local currents remain mesh/time sensitive. Ionization, preparation, physical edge and transport validity are unresolved. The displayed field peak includes traced exterior samples.")
    assets = Path(__file__).resolve().parents[1]/'visualizations/material_fields'
    html = (assets/'index.html').read_text()
    html = html.replace('This run starts with closed magnetic flux. It does not demonstrate FRC formation or particle confinement. External coils are not included; spatial accuracy and the plasma edge remain under study.',
        'A steady bias is followed by an opposing coil pulse. Colored rings show independently prescribed coil groups; grey lines show the tube. This coarse pre-ionized model does not establish FRC formation, physical transport or confinement.')
    if not d['success']:
        from html import escape
        html = html.replace('<div class="toolbar">', '<p style="color:#ffb697">'+escape(status)+'</p><div class="toolbar">')
    html = html.replace('<div class="toolbar">', '<p id="drive"></p><p id="drive-groups"></p><p id="motion"></p><div class="toolbar">')
    html = html.replace('material_vector_3d.png', 'material_experiment_3d.png')
    html = html.replace('Open the scientific image</a>', 'Open the scientific image</a> · <a href="material_experiment_sections.png">Signed fields and flux contours</a>')
    html = html.replace('/* DATA */', 'const DATA='+json.dumps(data, separators=(',', ':'))+';').replace('/* VIEWER */', (assets/'viewer.js').read_text())
    (out/'index.html').write_text(html); (out/'display_data.json').write_text(json.dumps(data, separators=(',', ':'))+'\n')
    gaps = np.diff([f['time_us'] for f in frames])*1000
    (out/'source.json').write_text(json.dumps(dict(source=str(source), inputs=inputs,
        display_interval_target_ns=frame_interval_ns, frame_count=len(frames), external_frame_files=external_frames,
        actual_gap_range_ns=[float(gaps.min()), float(gaps.max())] if len(gaps) else [],
        scope='Accepted saved states only; terminal status is shown explicitly. Full plasma and exterior B, interpolated field-line paths, sampled actual J. No temporal interpolation, topology certificate or particle trajectories.'), indent=2)+'\n')
    print('experiment display complete', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('source'); parser.add_argument('destination')
    parser.add_argument('--allow-partial', action='store_true')
    parser.add_argument('--frame-interval-ns', type=float)
    args = parser.parse_args(); build(args.source, args.destination, allow_partial=args.allow_partial, frame_interval_ns=args.frame_interval_ns)
