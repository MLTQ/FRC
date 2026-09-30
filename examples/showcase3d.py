"""Build a standalone showcase and optionally serve it with a narrow PNG export."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path


def build():
    root=Path(__file__).resolve().parents[1];assets=root/'visualizations/frc_showcase';output=root/'outputs/frc_showcase'
    html=(assets/'index.html').read_text().replace('/* DATA */','const DATA='+(output/'data.json').read_text()+';')
    html=html.replace('/* RENDERER */',(assets/'renderer.js').read_text()).replace('/* CONTROLLER */',(assets/'controller.js').read_text())
    (output/'index.html').write_text(html)
    return output


class Handler(SimpleHTTPRequestHandler):
    def do_POST(self):
        length=int(self.headers.get('Content-Length','0'))
        if self.path!='/save-poster' or not 0<length<30_000_000:
            self.send_error(400);return
        data=self.rfile.read(length)
        if not data.startswith(b'\x89PNG\r\n\x1a\n'):
            self.send_error(400);return
        (Path(self.directory)/'poster.png').write_bytes(data)
        self.send_response(200);self.end_headers();self.wfile.write(b'Saved')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--serve',action='store_true');parser.add_argument('--port',type=int,default=8765)
    args=parser.parse_args();output=build();print(output/'index.html',flush=True)
    if args.serve:
        ThreadingHTTPServer(('127.0.0.1',args.port),partial(Handler,directory=str(output))).serve_forever()
