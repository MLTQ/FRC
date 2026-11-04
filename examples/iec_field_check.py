"""End-to-end independent analytic electrostatic/particle refinement runs."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.constants import epsilon_0, elementary_charge as qe, electron_mass

from plasma_bench.iec_field import Field
from plasma_bench.iec_mesh import build_mesh


def run(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    results = dict(sphere=[], shell=[], oscillator=[])
    for geometry in ('sphere', 'shell'):
        for size in ((.04, .028, .02, .014, .01) if geometry == 'sphere' else (.04, .028, .02)):
            start = time.monotonic()
            mesh = build_mesh(out/f'{geometry}_{size}.npz', size=size, wire_size=size/2,
                              geometry=geometry)
            field = Field(mesh)
            radius = np.linalg.norm(field.nodes, axis=1)
            centres = field.nodes[field.tets].mean(axis=1)
            r = np.linalg.norm(centres, axis=1)
            if geometry == 'sphere':
                rho = 1e-9
                charge = np.bincount(field.tets.ravel(), np.repeat(rho*field.volume/4, 4), minlength=len(field.nodes))
                phi, electric, _, energy, _ = field.solve(charge, 0.)
                analytic_phi = rho*(.2**2-radius**2)/(6*epsilon_0)
                analytic_e = rho*centres/(3*epsilon_0)
                analytic_energy = 2*np.pi*rho**2*.2**5/(45*epsilon_0)
            else:
                phi, electric, _, energy, _ = field.solve(np.zeros(len(field.nodes)), 1.)
                a, b = .075, .2
                analytic_phi = (1/radius-1/b)/(1/a-1/b)
                analytic_e = centres/(r[:, None]**3*(1/a-1/b))
                analytic_energy = .5*4*np.pi*epsilon_0*a*b/(b-a)
            error_e = np.sqrt(np.sum(field.volume[:, None]*(electric-analytic_e)**2)/
                              np.sum(field.volume[:, None]*analytic_e**2))
            row = dict(size_m=size, nodes=len(field.nodes), tetrahedra=len(field.tets),
                       potential_relative_l2=float(np.linalg.norm(phi-analytic_phi)/np.linalg.norm(analytic_phi)),
                       field_volume_relative_l2=float(error_e), field_energy_relative_error=float(abs(energy/analytic_energy-1)),
                       elapsed_s=time.monotonic()-start)
            results[geometry].append(row)
            print(geometry, row, flush=True)
            if geometry == 'sphere':
                # Independent uniform positive charge has an exact electron oscillator.
                omega = np.sqrt(qe*rho/(3*epsilon_0*electron_mass))
                x0 = np.array([[.031, .009, .012]])
                for subdivisions in (100, 200):
                    dt = 2*np.pi/omega/subdivisions
                    x, v = x0.copy(), np.zeros_like(x0)
                    cells = field.locate(x)
                    k0 = .5*electron_mass*omega**2*np.sum(x0*x0)
                    max_energy = 0.
                    for _ in range(subdivisions):
                        v -= .5*dt*qe/electron_mass*electric[cells]
                        x, cells, hit = field.drift(x, x+dt*v, cells)
                        if hit.any():
                            raise RuntimeError('Oscillator left the gas unexpectedly')
                        v -= .5*dt*qe/electron_mass*electric[cells]
                        h = .5*electron_mass*(np.sum(v*v)+omega**2*np.sum(x*x))
                        max_energy = max(max_energy, abs(h/k0-1))
                    results['oscillator'].append(dict(size_m=size, steps=subdivisions,
                        position_relative_error=float(np.linalg.norm(x-x0)/np.linalg.norm(x0)),
                        maximum_analytic_energy_relative_error=max_energy))
            (out/'report.json').write_text(json.dumps(results, indent=2, allow_nan=False))
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    run(parser.parse_args().output)
