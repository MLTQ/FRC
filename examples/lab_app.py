"""Loopback-only experiment app; owns only the IEC child runs it launches."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import secrets
import signal
import subprocess
import sys
import threading
import time
from urllib.parse import parse_qs, unquote, urlparse

from plasma_bench.lab_catalog import catalog, runs, validate
from plasma_bench.iec_observables import ObservableStore


class Lab:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.token = secrets.token_urlsafe(24)
        self.lock = threading.Lock()
        self.children = {}
        self.observables = ObservableStore(self.root)
        (self.root/'outputs/lab_runs').mkdir(parents=True, exist_ok=True)

    def launch(self, values):
        config = validate(values)
        with self.lock:
            if any(child.poll() is None for child in self.children.values()):
                raise ValueError('A local electrostatic run is already active. Wait or stop it before starting another.')
            identifier = time.strftime('iec_%Y%m%d_%H%M%S')+'_'+secrets.token_hex(2)
            base = self.root/'outputs/lab_runs'
            config_path = base/(identifier+'.json')
            config_path.write_text(json.dumps(config, allow_nan=False))
            log = (base/(identifier+'.log')).open('w')
            environment = os.environ | dict(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
            command = [sys.executable, '-m', 'examples.iec_3d', str(base/identifier), '--config', str(config_path)]
            # Reuse only a cache with exactly matching declared geometry/resolution.
            cache = self.root/'outputs/iec_mesh_dev/cage.npz'
            if config['mesh_mm'] == 18 and config['wire_mesh_mm'] == 1 and cache.exists():
                command += ['--mesh', str(cache)]
            child = subprocess.Popen(command, cwd=self.root, env=environment, stdout=log, stderr=subprocess.STDOUT)
            log.close()
            self.children[identifier] = child
            return dict(id=identifier, status='starting')

    def cancel(self, identifier):
        with self.lock:
            child = self.children.get(identifier)
            if child is None or child.poll() is not None:
                raise ValueError('This app has no active owned run with that identifier')
            directory = self.root/'outputs/lab_runs'/identifier
            if not directory.exists():
                raise ValueError('Run is still starting; try again shortly')
            (directory/'STOP').touch()
            return dict(id=identifier, status='stop requested')


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, lab=None, **kwargs):
        self.lab = lab
        super().__init__(*args, directory=str(lab.root/'visualizations/lab'), **kwargs)

    def log_message(self, fmt, *args):
        if args and str(args[0]).startswith('POST'):
            super().log_message(fmt, *args)

    def local_host(self):
        return self.headers.get('Host', '').lower() in {
            f'localhost:{self.server.server_port}', f'127.0.0.1:{self.server.server_port}'}

    def json(self, value, status=200):
        body = json.dumps(value, allow_nan=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if not self.local_host():
            return self.send_error(403, 'Loopback host required')
        path = urlparse(self.path).path
        if path == '/api/state':
            return self.json(dict(experiments=catalog(), runs=runs(self.lab.root), token=self.lab.token))
        if path == '/api/observables':
            try:
                query = parse_qs(urlparse(self.path).query)
                identifier = query.get('run', [''])[0]
                index = int(query.get('frame', ['-1'])[0])
                magnetic = query.get('magnetic', ['0'])[0]
                if magnetic not in ('0', '1'):
                    raise ValueError('Magnetic diagnostic must be zero or one')
                return self.json(self.lab.observables.get(identifier, index, magnetic == '1'))
            except (ValueError, OSError, KeyError) as error:
                return self.json(dict(error=str(error)), 400)
        if path.startswith('/data/'):
            relative = unquote(path[len('/data/'):])
            # Only scientific output files; never repo secrets/source or directory lists.
            file = (self.lab.root/'outputs'/relative).resolve()
            if not file.is_relative_to(self.lab.root/'outputs') or not file.is_file():
                return self.send_error(404)
            if file.suffix.lower() not in ('.json', '.html', '.js', '.css', '.bin', '.f16', '.f32', '.rgba16f', '.png', '.jpg', '.svg'):
                return self.send_error(404)
            self.send_response(200)
            self.send_header('Content-Type', self.guess_type(str(file)))
            self.send_header('Content-Length', str(file.stat().st_size))
            self.end_headers()
            try:
                with file.open('rb') as stream:
                    import shutil
                    shutil.copyfileobj(stream, self.wfile)
            except (BrokenPipeError, ConnectionResetError):
                pass
            return
        if path.startswith('/api/'):
            return self.send_error(404)
        if path not in ('/', '/index.html', '/app.js', '/scene.js', '/fields.js', '/style.css'):
            return self.send_error(404)
        super().do_GET()

    def do_POST(self):
        # A browser on another site cannot launch local simulations cross-origin.
        origin = self.headers.get('Origin')
        if not self.local_host() or self.headers.get('X-Lab-Token') != self.lab.token or (origin and origin != 'http://'+self.headers.get('Host', '')):
            return self.json(dict(error='Local session token required'), 403)
        try:
            length = int(self.headers.get('Content-Length', '0'))
            if not 0 < length <= 16000:
                raise ValueError('Invalid request size')
            value = json.loads(self.rfile.read(length))
            path = urlparse(self.path).path
            if path == '/api/launch':
                return self.json(self.lab.launch(value))
            if path == '/api/stop':
                return self.json(self.lab.cancel(value.get('id')))
            return self.send_error(404)
        except (ValueError, TypeError, KeyError) as error:
            self.json(dict(error=str(error)), 400)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8792)
    args = parser.parse_args()
    lab = Lab(Path(__file__).resolve().parents[1])
    server = ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, lab=lab))
    print(f'Experiment lab: http://localhost:{args.port}', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        # Running physical children retain their output and finish independently.


if __name__ == '__main__':
    main()
