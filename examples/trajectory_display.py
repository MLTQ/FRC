"""Assemble an offline display of the current archived simulation state."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from .implicit_plasma_records import sha


def build():
    root = Path(__file__).resolve().parents[1]
    source = root/'outputs/trajectory_display_data'
    output = root/'outputs/trajectory_display'; output.mkdir(exist_ok=True)
    summary = json.loads((source/'summary.json').read_text())
    if not summary['terminal'] or sha(source/'data.json') != summary['data_sha256']:
        raise ValueError('Complete unchanged presentation data required')
    assets = root/'visualizations/trajectory_display'
    html = (assets/'index.html').read_text()
    for marker, path, prefix, suffix in [('DATA', source/'data.json', 'const DATA=', ';'),
                                        ('RENDERER', assets/'renderer.js', '', ''),
                                        ('CONTROLLER', assets/'controller.js', '', '')]:
        html = html.replace('/* '+marker+' */', prefix+path.read_text()+suffix)
    target = output/'index.html'; target.write_text(html)
    names = [*sorted(assets.glob('*')), Path(__file__), Path(__file__).with_suffix('.md'), source/'data.json', source/'summary.json']
    manifest = dict(inputs={str(p.relative_to(root)): sha(p) for p in names},
                    html_sha256=sha(target), source_trajectory=summary['scope'])
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return output


class Handler(SimpleHTTPRequestHandler):
    def do_POST(self):
        try: length = int(self.headers.get('Content-Length', '0'))
        except ValueError: self.send_error(400); return
        if self.path != '/save-poster' or not 0 < length < 30_000_000:
            self.send_error(400); return
        data = self.rfile.read(length)
        if len(data) != length or not data.startswith(b'\x89PNG\r\n\x1a\n'):
            self.send_error(400); return
        (Path(self.directory)/'poster.png').write_bytes(data)
        self.send_response(200); self.end_headers(); self.wfile.write(b'Saved')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--serve', action='store_true')
    parser.add_argument('--port', type=int, default=8786); args = parser.parse_args()
    directory = build(); print(directory/'index.html', flush=True)
    if args.serve:
        ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, directory=str(directory))).serve_forever()
