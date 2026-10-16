"""Dense, read-only current volumes for an existing saved-state animation."""
import argparse
import json
import os
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from plasma_bench.material_experiment import build as experiment
from plasma_bench.material_plasma_samples import material_samples
from .material_coil_plot import state_at
from .implicit_plasma_records import sha


def volume(model, state, plasma, metadata, time, extent):
    r, z = np.meshgrid(np.linspace(1e-5, 1-1e-5, 65)*plasma['radius_m'],
                       np.linspace(-1+1e-5, 1-1e-5, 257)*plasma['half_length_m'], indexing='ij')
    density = metadata['central_density_m3']*((1-(r/plasma['radius_m'])**2)*
                                             (1-(z/plasma['half_length_m'])**2))**1.5
    response = model.response(state, rate=False, time=time)
    # J depends on the full response, but has no additive prescribed-source term.
    # Avoid computing unused coil B at every volume point; the reconstructed
    # drift still comes from the unchanged complete physical response.
    interior = SimpleNamespace(space=model.space, fluid=model.fluid)
    samples = material_samples(interior, state, r, z, density, time=time, response=response)
    points = samples['points'].reshape(-1, 2)
    currents = samples['current'].reshape(-1, 3)
    # Strictly interior interpolation; no extrapolated current or artificial halo.
    rmax = points[:, 0].max(); zmin = points[:, 1].min(); zmax = points[:, 1].max()
    rr, zz = np.meshgrid(np.linspace(0, rmax, 96), np.linspace(zmin, zmax, 384), indexing='ij')
    values = LinearNDInterpolator(points, currents)(rr, zz)
    supported = np.isfinite(values).all(axis=-1)
    # A convex interpolation hull can span a concave material boundary. Explicitly
    # keep only points inside the actual polygonal Q1 domain as well.
    from matplotlib.path import Path as Polygon
    coords = state.fluid.coordinates[model.fluid.domain.ids]
    boundary = np.concatenate((coords[:, 0], coords[-1, 1:], coords[-2::-1, -1], coords[0, -2:0:-1]))
    supported &= Polygon(boundary).contains_points(np.column_stack((rr.ravel(), zz.ravel()))).reshape(rr.shape)
    texture = np.zeros(rr.shape+(4,), dtype='<f4')
    texture[..., :3] = np.where(supported[..., None], values, 0)
    texture[..., 3] = supported
    maxima = np.r_[np.linalg.norm(currents, axis=1).max(), abs(currents).max(axis=0)]
    return texture, dict(shape=[384, 96], bounds=[0., rmax/extent, zmin/extent, zmax/extent],
                         jmax=maxima.tolist(), supported_texels=int(supported.sum()))


def build(existing, destination, limit=None):
    original = Path(existing).resolve(); out = Path(destination).resolve()
    data = json.loads((original/'display_data.json').read_text())
    provenance = json.loads((original/'source.json').read_text())
    source = Path(provenance['source']).resolve()
    summary = json.loads((source/'summary.json').read_text())
    if out.exists() or not summary['terminal'] or not summary['success']:
        raise ValueError('Fresh output and a completed source required')
    fluid, model, metadata = experiment(summary['config'])
    out.mkdir(parents=True); (out/'volumes').mkdir()
    (out/'.gitignore').write_text('volumes/\nframes\n')
    (out/'frames').symlink_to(os.path.relpath(original/'frames', out), target_is_directory=True)
    for filename in ('material_experiment_3d.png', 'material_experiment_sections.png'):
        (out/filename).symlink_to(os.path.relpath(original/filename, out))
    frames = data['frames'][:limit] if limit else data['frames']
    for i, frame in enumerate(frames):
        state = state_at(source, frame['step'])
        texture, description = volume(model, state, summary['config']['plasma'], metadata,
                                      frame['step']*summary['dt'], data['physical_scale_m'])
        filename = f"volumes/step_{frame['step']:05}.f32"
        texture.tofile(out/filename)
        description.update(url=filename, sha256=sha(out/filename))
        frame['volume'] = description
        frame['jmax'] = np.maximum(frame['jmax'], description['jmax']).tolist()
        print(f"volume {i+1}/{len(frames)} step {frame['step']}", flush=True)
    data.update(frames=frames, detailed=True, initial_yaw=-1.14, initial_pitch=.3,
                initial_zoom=4.4, initial_frame=int(np.argmax([f['jmax'][0] for f in frames])),
                jmax=np.max([f['jmax'] for f in frames], axis=0).tolist())
    data['scope'] += (' Detailed current volume samples the actual material basis at 65×257 labels, '
        'then interpolates onto a 96×384 meridional display texture. Empty exterior is masked. '
        'Field-line density uses exact azimuthal symmetry of the stored paths. '
        'Rendering resolution does not increase the 4×16 simulation resolution. '
        'Current maxima include the dense basis samples; opacity and threshold affect display only.')
    assets = Path(__file__).resolve().parents[1]/'visualizations/material_fields'
    html = (assets/'volume.html').read_text().replace('/* DATA */', 'const DATA='+json.dumps(data, separators=(',', ':'))+';')
    html = html.replace('/* VOLUME */', (assets/'volume.js').read_text()).replace('/* VIEWER */', (assets/'viewer.js').read_text())
    (out/'index.html').write_text(html)
    (out/'display_data.json').write_text(json.dumps(data, separators=(',', ':'))+'\n')
    provenance.update(original_display=str(original), original_display_sha256=sha(original/'display_data.json'),
                      volume_frames=len(frames), complete_display=limit is None,
                      volume_method='Compatible material samples; linear physical interpolation; axisymmetric GPU volume. No solver evolution.')
    (out/'source.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Detailed display complete:', out, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('existing'); parser.add_argument('destination'); parser.add_argument('--limit', type=int)
    args = parser.parse_args(); build(args.existing, args.destination, args.limit)
