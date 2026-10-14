"""Read-only presentation of archived B and signed current components."""
import argparse
import base64
from functools import partial
from http.server import ThreadingHTTPServer
import json
from pathlib import Path
import numpy as np

from .current_state_display_v2 import replace_once, verify_inputs
from .current_state_poster_v3 import render
from .implicit_plasma_records import sha
from .trajectory_display import Handler

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'outputs/trajectory_display_data'


def research_inputs():
    """Check the completed evidence actually cited by the status fragment."""
    checked = {}
    checkpoint_path = ROOT / 'outputs/species_budget_diagnostics/checkpoint_v2.json'
    checkpoint = json.loads(checkpoint_path.read_text())
    if not checkpoint['terminal']:
        raise ValueError('Complete species checkpoint required')
    checked[checkpoint_path] = sha(checkpoint_path)
    checked.update({ROOT / p: h for p, h in checkpoint['evidence_files'].items()
                    if p.startswith('outputs/')})
    provenance = ROOT / 'outputs/species_budget_diagnostics/active_source_v2/provenance.json'
    checked[provenance] = sha(provenance)
    checked.update({ROOT / p: h for p, h in json.loads(provenance.read_text())['source_files'].items()})
    for name in ('coupled_fluid_probe_F16', 'coupled_fluid_probe_F32', 'coupled_fluid_audit_F16'):
        verifier_path = ROOT / f'outputs/{name}_replay_verification/summary.json'
        record = json.loads(verifier_path.read_text())
        if not record['terminal'] or not record['bitwise_metrics']:
            raise ValueError('Complete scientific replay required')
        checked[verifier_path] = sha(verifier_path)
        for variant in ('original', 'repeated'):
            directory = ROOT / record[variant + '_directory']
            checked[directory / 'summary.json'] = record[variant + '_summary_sha256']
            for filename, item in record['arrays'].items():
                if not item['bitwise_equal']:
                    raise ValueError('Scientific raw replay mismatch')
                checked[directory / filename] = item[variant + '_sha256']
            prov = directory / 'provenance.json'
            checked[prov] = sha(prov)
            checked.update({ROOT / p: h for p, h in json.loads(prov.read_text())['source_files'].items()})
    for path, expected in checked.items():
        if sha(path) != expected:
            raise ValueError(f'Changed archived evidence: {path}')
    return checked


def package_currents(data, summary):
    """Pack |J| and the signed SI components; leave every magnetic path intact."""
    records = []
    for frame, measured in zip(data['frames'], summary['frames'], strict=True):
        with np.load(SOURCE / measured['arrays_file']) as raw:
            current = raw['J']
            magnitude = np.linalg.norm(current, axis=-1)
            fields = np.concatenate((magnitude[..., None], current), axis=-1)
        texture = fields.transpose(1, 0, 2).astype('<f4')
        old_texture = np.frombuffer(base64.b64decode(frame['texture']), dtype='<f4').reshape(texture.shape)
        if not np.array_equal(texture[..., 0], old_texture[..., 0]):
            raise ValueError('Original magnitude texture must remain bitwise identical')
        packed = base64.b64encode(texture.tobytes()).decode()
        decoded = np.frombuffer(base64.b64decode(packed), dtype='<f4').reshape(texture.shape)
        if not np.array_equal(decoded, texture):
            raise ValueError('Packed component round trip changed values')
        error = decoded.astype(float) - fields.transpose(1, 0, 2)
        absolute = np.max(np.abs(error), axis=(0, 1))
        limits = np.max(np.abs(fields), axis=(0, 1))
        relative = absolute / np.maximum(limits, np.finfo(float).tiny)
        if np.max(relative) > np.finfo(np.float32).eps:
            raise ValueError('Unexpected float32 packing error')
        frame['texture'] = packed
        frame['component_ranges'] = np.stack((fields.min(axis=(0, 1)), fields.max(axis=(0, 1))), axis=1).tolist()
        frame['component_absmax'] = limits.tolist()
        records.append(dict(step=frame['step'], component_ranges_A_m2=frame['component_ranges'],
                            max_float32_absolute_error_A_m2=absolute.tolist(),
                            max_float32_error_over_component_peak=relative.tolist(),
                            original_magnitude_texture_identical=True))
    data['component_absmax'] = np.max([f['component_absmax'] for f in data['frames']], axis=0).tolist()
    return records


def build(output):
    if output.exists():
        raise FileExistsError(f'Preserve completed presentation: {output}')
    summary, _, old_checked = verify_inputs()
    checked = {p: h for p, h in old_checked.items() if str(p.relative_to(ROOT)).startswith('outputs/')}
    checked.update(research_inputs())
    assets = ROOT / 'visualizations/trajectory_display'
    old_additions = ROOT / 'visualizations/current_state_display'
    additions = ROOT / 'visualizations/current_state_display_v3'
    direct = [*sorted(additions.glob('*')), assets / 'index.html', old_additions / 'status.css',
              old_additions / 'controls.js', SOURCE / 'summary.json', Path(__file__),
              Path(__file__).with_suffix('.md'), ROOT / 'examples/current_state_poster_v3.py',
              ROOT / 'examples/current_state_poster_v3.md', ROOT / 'examples/current_state_display_v2.py',
              ROOT / 'examples/trajectory_display.py', ROOT / 'examples/implicit_plasma_records.py']
    checked.update({p: sha(p) for p in direct})
    data = json.loads((SOURCE / 'data.json').read_text())
    packing = package_currents(data, summary)
    html = (assets / 'index.html').read_text()
    replacements = [
        ('FRC sandbox · Current simulation fields', 'FRC sandbox · Magnetic fields and signed currents'),
        ('The current simulation fields.', 'Magnetic fields. Plasma currents.'),
        ('Current density and magnetic field lines from the latest evolved state.',
         'Explore the latest accepted state and three earlier saved times.'),
        ('</style>', (old_additions / 'status.css').read_text() +
         'select{background:#10232f;color:#d5e8e3;border:1px solid #47616c;border-radius:6px;padding:8px}'
         '.metric b{font-size:clamp(14px,1.6vw,21px)}.section-image{display:block;width:100%;margin:16px 0}'
         '</style>'),
        ('<p>Volume color = |J|', '<p id="current-description">Volume color = |J|'),
        ('<div class="bar warm">', '<div id="current-bar" class="bar warm">'),
        ('<span class="small">SAMPLED PEAK |J|', '<span id="current-metric-label" class="small">SAMPLED PEAK |J|'),
        ('<button id="reset">Reset view</button>', '<button id="reset">Reset view</button>'
         '<div class="group" aria-label="Camera presets"><button data-view="perspective">3D</button>'
         '<button data-view="side">Side</button><button data-view="end">End</button></div>'),
        ('<button id="scale"', '<label class="hint" for="component">Current view</label>'
         '<select id="component"><option value="0">Magnitude |J|</option><option value="1">Radial Jr</option>'
         '<option value="2">Azimuthal Jθ</option><option value="3">Axial Jz</option></select><button id="scale"'),
        ('<div class="readouts">', '<p class="hint" id="component-note"></p>'
         '<div class="readouts">'),
        ('<div class="footer">', (additions / 'status.html').read_text() +
         '<details><summary>Scientific image and signed cross-sections</summary>'
         '<a href="current_state_3d.png" target="_blank"><img class="section-image" src="current_state_3d.png" '
         'alt="Last accepted state at 28.33 microseconds: 3D current magnitude and magnetic field lines, '
         'signed axial current and axial field departure cross-sections"></a></details><div class="footer">'),
        ('<span id="status" role="status">', '<a href="current_state_3d.png" target="_blank">Open scientific image ↗</a>'
         '<span id="status" role="status">')]
    for old, new in replacements:
        html = replace_once(html, old, new)
    controller = (additions / 'controller.js').read_text() + '\n' + (old_additions / 'controls.js').read_text()
    for marker, content in [('DATA', 'const DATA=' + json.dumps(data, separators=(',', ':')) + ';'),
                            ('RENDERER', (additions / 'renderer.js').read_text()), ('CONTROLLER', controller)]:
        html = replace_once(html, '/* ' + marker + ' */', content)
    output.mkdir(parents=True)
    (output / 'index.html').write_text(html)
    render(SOURCE, summary, data, output / 'current_state_3d.png')
    (output / 'packing.json').write_text(json.dumps(dict(frames=packing), indent=2) + '\n')
    for path, expected in checked.items():
        if sha(path) != expected:
            raise ValueError(f'Input changed during presentation build: {path}')
    manifest = dict(terminal=True, latest_accepted_step=17, latest_accepted_time_us=summary['sampled_time_max_us'],
                    simulation_advanced=False, unit_tests_run=False, scope=summary['scope'],
                    inputs={str(p.relative_to(ROOT)): h for p, h in sorted(checked.items())},
                    outputs={p.name: sha(p) for p in sorted(output.iterdir())},
                    browser_verification='Unavailable: Mac locked. Static rendering does not verify WebGL interaction.')
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/current_state_display_v3')
    parser.add_argument('--serve', action='store_true')
    parser.add_argument('--port', type=int, default=8788)
    args = parser.parse_args()
    output = build(args.output.resolve())
    print(output / 'index.html', flush=True)
    if args.serve:
        ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, directory=str(output))).serve_forever()


if __name__ == '__main__':
    main()
