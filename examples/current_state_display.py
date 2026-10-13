"""Repackage the archived 3D fields with current research status, read-only."""
import argparse
from functools import partial
from http.server import ThreadingHTTPServer
import json
from pathlib import Path

from .implicit_plasma_records import sha
from .trajectory_display import Handler
from .current_state_poster import render

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / 'outputs/trajectory_display_data'
AUDIT_DIR = ROOT / 'outputs/trajectory_spatial_parts_diagnostics'


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected exactly one frozen presentation marker: {old[:80]}')
    return text.replace(old, new, 1)


def verify_inputs():
    summary = json.loads((DATA_DIR / 'summary.json').read_text())
    checkpoint = json.loads((AUDIT_DIR / 'checkpoint.json').read_text())
    if not summary['terminal'] or not checkpoint['terminal']:
        raise ValueError('Complete display data and audits required')
    checked = {DATA_DIR / 'data.json': summary['data_sha256']}
    for frame in summary['frames']:
        checked[DATA_DIR / frame['arrays_file']] = frame['arrays_sha256']
        checked[ROOT / frame['source_file']] = frame['source_sha256']
    for name, digest in checkpoint['evidence_files'].items():
        checked[ROOT / name] = digest
    # Root status/docs may advance after the historical research checkpoint;
    # its scientific outputs must retain every recorded byte.
    for path, digest in checked.items():
        if str(path.relative_to(ROOT)).startswith('outputs/') and sha(path) != digest:
            raise ValueError(f'Archived evidence changed: {path}')
    return summary, checkpoint, checked


def build(output):
    if output.exists():
        raise FileExistsError(f'Preserve completed output: {output}')
    summary, checkpoint, checked = verify_inputs()
    assets = ROOT / 'visualizations/trajectory_display'
    additions = ROOT / 'visualizations/current_state_display'
    html = (assets / 'index.html').read_text()
    html = replace_once(html, 'FRC sandbox · Current simulation fields', 'FRC sandbox · Latest fields and audit status')
    html = replace_once(html, 'The current simulation fields.', 'The latest saved plasma fields.')
    html = replace_once(html, 'Current density and magnetic field lines from the latest evolved state.',
                        'Explore the last accepted state, then compare the four saved times.')
    html = replace_once(html, '</style>', (additions / 'status.css').read_text() + '</style>')
    html = replace_once(html, '<div class="readouts">',
                        (additions / 'status.html').read_text() + '<div class="readouts">')
    html = replace_once(html, '<button id="reset">Reset view</button>',
                        '<button id="reset">Reset view</button><div class="group" aria-label="Camera presets">'
                        '<button data-view="perspective">3D</button><button data-view="side">Side</button>'
                        '<button data-view="end">End</button></div>')
    html = replace_once(html, '<span id="status" role="status">',
                        '<a href="current_state_3d.png" target="_blank">Open scientific image ↗</a>'
                        '<span id="status" role="status">')
    controller = (assets / 'controller.js').read_text().replace(
        'Saved poster.png in outputs/trajectory_display', 'Saved poster.png beside this viewer')
    controller += '\n' + (additions / 'controls.js').read_text()
    for marker, value in [('DATA', 'const DATA=' + (DATA_DIR / 'data.json').read_text() + ';'),
                          ('RENDERER', (assets / 'renderer.js').read_text()), ('CONTROLLER', controller)]:
        html = replace_once(html, '/* ' + marker + ' */', value)
    output.mkdir(parents=True)
    (output / 'index.html').write_text(html)
    data = json.loads((DATA_DIR / 'data.json').read_text())
    render(DATA_DIR, summary, data, output / 'current_state_3d.png')
    inputs = [*sorted(assets.glob('*')), *sorted(additions.glob('*')),
              Path(__file__), Path(__file__).with_suffix('.md'),
              ROOT / 'examples/current_state_poster.py', ROOT / 'examples/current_state_poster.md',
              ROOT / 'examples/trajectory_display.py', DATA_DIR / 'data.json', DATA_DIR / 'summary.json',
              AUDIT_DIR / 'checkpoint.json', *[DATA_DIR / f['arrays_file'] for f in summary['frames']]]
    manifest = dict(terminal=True, scope=summary['scope'],
                    latest_accepted_step=summary['frames'][-1]['step'],
                    latest_accepted_time_us=summary['sampled_time_max_us'],
                    simulation_advanced=False, audit_source_sha256=checkpoint['source_sha256'],
                    verified_output_inputs=sum(str(p.relative_to(ROOT)).startswith('outputs/') for p in checked),
                    inputs={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                    outputs={p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()},
                    browser_visual_check='Unavailable: Mac locked; no interactive verification claimed.',
                    unit_tests_run=False)
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/current_state_display')
    parser.add_argument('--serve', action='store_true')
    parser.add_argument('--port', type=int, default=8787)
    args = parser.parse_args()
    output = build(args.output.resolve())
    print(output / 'index.html', flush=True)
    if args.serve:
        ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, directory=str(output))).serve_forever()


if __name__ == '__main__':
    main()
