"""Independent electrostatic comparison and full open-cage vacuum calculation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0

from plasma_bench.iec_vacuum import (Sources, cage_discretization, directions,
                                     fit_conductor, solve_cage, truncated_cage)


def errors(solution, points, weights=None):
    error = solution.evaluate(points)/solution.voltage-1
    return dict(relative_rms=float(np.sqrt(np.average(error**2, weights=weights))),
                relative_max=float(np.max(np.abs(error))))


def spheres():
    """Compare computed fields with the independent spherical Gauss-law solution."""
    a, wall, voltage = .075, .2, 1.
    unit = directions(127, .71)
    radii = np.array([1.01*a, .5*(a+wall), .99*wall])
    x = np.concatenate([r*unit for r in radii])
    r = np.linalg.norm(x, axis=1)
    exact_phi = voltage*(1/r-1/wall)/(1/a-1/wall)
    exact_e = voltage*x/(r[:, None]**3*(1/a-1/wall))
    exact_c = 4*np.pi*epsilon_0/(1/a-1/wall)
    rows = []
    for n in (12, 24, 48, 96):
        source = Sources(np.empty((0, 3)), np.empty((0, 3)), .65*a*directions(n), wall)
        solution = fit_conductor(source, a*directions(3*n, .29), voltage=voltage)
        phi, field = solution.evaluate(x, field=True)
        row = dict(source_points=n, potential_error_max_over_voltage=float(np.max(abs(phi-exact_phi))/voltage),
                   field_relative_l2=float(np.linalg.norm(field-exact_e)/np.linalg.norm(exact_e)),
                   capacitance_relative_error=float(abs(solution.capacitance/exact_c-1)),
                   off_fit_surface=errors(solution, a*directions(5*n+1, .83)), fit=solution.fit)
        rows.append(row)
        print('sphere', json.dumps(row), flush=True)
    return dict(assumptions='Concentric spherical conductors; reference is independent Gauss law.',
                exact_capacitance_f=exact_c, rows=rows)


def finite_line():
    """Independent bisector/axial field of a finite uniform line in free space."""
    length, wall = .02, 1e6
    sources = Sources(np.array([[0., 0., -length/2]]), np.array([[0., 0., length/2]]),
                      np.empty((0, 3)), wall, quadrature=9)
    r = np.array([.001, .003, .01, .03])
    x = np.column_stack((r, np.zeros((len(r), 2))))
    phi, electric = sources.kernel(x, field=True)
    exact_phi = 2/length*np.log((np.sqrt(r*r+(length/2)**2)+length/2)/r)-1/wall
    exact_e = 1/(r*np.sqrt(r*r+(length/2)**2))
    axial = np.array([[0., 0., .03], [0., 0., -.03]])
    _, axis_e = sources.kernel(axial, field=True)
    axis_exact = np.sign(axial[:, 2])/(axial[:, 2]**2-(length/2)**2)
    return dict(potential_relative_max=float(np.max(abs(phi[:, 0]/exact_phi-1))),
                bisector_field_relative_max=float(np.max(abs(electric[:, 0, 0]/exact_e-1))),
                axial_field_relative_max=float(np.max(abs(axis_e[:, 0, 2]/axis_exact-1))),
                assumption='Free line analytic reference; finite grounded wall correction beyond constant is negligible at R=1e6m.')


def interior_field(solution, radius):
    rays = directions(61, .17)
    radial = []
    for fraction in (.25, .5, .75):
        r = fraction*radius
        phi, field = solution.evaluate(r*rays, field=True)
        er = np.sum(field*rays, axis=1)/solution.voltage
        radial.append(dict(radius_m=r, minimum_outward_E_per_bias=float(er.min()),
                           maximum_outward_E_per_bias=float(er.max()),
                           angular_potential_span_over_bias=float(np.ptp(phi)/solution.voltage)))
    point = radius*np.array([.25, .13, .07])
    divergence = []
    for h in (.02*radius, .01*radius):
        offsets = h*np.eye(3)
        _, field = solution.evaluate(np.vstack((point+offsets, point-offsets)), field=True)
        div = np.trace(field[:3]-field[3:])/(2*h)
        divergence.append(dict(spacing_m=h, divergence_times_radius_squared_over_bias=float(div*radius**2/solution.voltage)))
    return dict(radial_samples=radial, vacuum_divergence_refinement=divergence,
                interpretation='Both signs of radial E mean no inward force in every direction for positive ions. Centre E=0 alone would not establish or exclude a trap.')


def plot(cage, solution, out):
    fig = plt.figure(figsize=(14, 4.5), layout='constrained')
    ax = fig.add_subplot(131, projection='3d')
    for ia, ib in cage.edges:
        v = 100*cage.vertices[[ia, ib]]
        ax.plot(*v.T, color='#b64b2e', lw=1.5)
    ax.set(xlabel='x (cm)', ylabel='y (cm)', zlabel='z (cm)',
           title='Open wire cage · all 90 edges')
    ax.set_box_aspect((1, 1, 1))
    ax = fig.add_subplot(132)
    axis = np.linspace(-.11, .11, 151)
    xx, zz = np.meshgrid(axis, axis)
    points = np.column_stack((xx.ravel(), np.zeros(xx.size), zz.ravel()))
    outside = cage.wire_distance(points) > cage.wire_radius
    phi = np.full(len(points), np.nan)
    phi[outside] = solution.evaluate(points[outside])/solution.voltage
    im = ax.pcolormesh(100*xx, 100*zz, phi.reshape(xx.shape), shading='nearest', cmap='viridis')
    fig.colorbar(im, ax=ax, label='Potential / cage voltage')
    ax.set(xlabel='x (cm)', ylabel='z (cm)', aspect='equal', title='Vacuum potential · y = 0')
    ax = fig.add_subplot(133)
    radius = np.linspace(0, .19, 301)
    vectors = np.vstack((np.eye(3), np.ones(3)/np.sqrt(3), cage.cap_directions[:2]))
    for k, vector in enumerate(vectors):
        points = radius[:, None]*vector
        outside = cage.wire_distance(points) > cage.wire_radius
        phi = np.full(len(points), np.nan)
        phi[outside] = solution.evaluate(points[outside])/solution.voltage
        ax.plot(100*radius, phi, label=('x', 'y', 'z', 'diagonal', 'ray 5', 'ray 6')[k])
    ax.set(xlabel='Distance from centre (cm)', ylabel='Potential / cage voltage',
           title='Separate 3D rays · no angular average')
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    cage_cm = 200*(np.linalg.norm(cage.vertices[0])+cage.wire_radius)
    fig.suptitle(f'IEC reference · cage {cage_cm:g} cm, wire {2000*cage.wire_radius:g} mm · grounded spherical chamber {200*cage.wall_radius:g} cm\n'
                 'Vacuum field only: no plasma, prescribed trap, heating or fusion result', fontsize=12)
    fig.savefig(out/'vacuum_field.png', dpi=180)
    fig.savefig(out/'vacuum_field.pdf')
    plt.close(fig)


def main(out, panels=(2, 4, 6), diameter=.15, wire=.002, wall=.2):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    saved = out/'source'
    saved.mkdir()
    project = Path(__file__).resolve().parents[1]
    hashes = {}
    for folder in ('plasma_bench', 'examples'):
        for suffix in ('py', 'md'):
            source = project/folder/f'iec_vacuum.{suffix}'
            shutil.copyfile(source, saved/f'{folder}_iec_vacuum.{suffix}')
            hashes[str(source.relative_to(project))] = hashlib.sha256(source.read_bytes()).hexdigest()
    (out/'source_hashes.json').write_text(json.dumps(hashes, indent=1))
    inputs = dict(cage_outside_diameter_m=diameter, wire_diameter_m=wire, wall_radius_m=wall,
                  reference_boundary='grounded sphere, no support/feedthrough', plasma='none', panels=list(panels))
    (out/'inputs.json').write_text(json.dumps(inputs, indent=1))
    reference = spheres()
    line = finite_line()
    print('finite_line', json.dumps(line), flush=True)
    (out/'spheres.json').write_text(json.dumps(reference, indent=1))
    cage = truncated_cage(diameter, wire, wall)
    _, holdout, weights = cage_discretization(cage, max(panels), axial_shift=.23, holdout=True)
    rows = []
    for n in panels:
        solution = solve_cage(cage, n)
        centre_phi, centre_e = solution.evaluate([[0., 0., 0.]], field=True)
        wall_phi = solution.evaluate(wall*directions(131, .39))
        row = dict(panels_per_edge=n, fit=solution.fit, off_fit_surface=errors(solution, holdout, weights),
                   centre_potential_over_bias=float(centre_phi[0]),
                   centre_field_per_bias_V_m_per_V=centre_e[0].tolist(),
                   wall_potential_error_over_bias=float(np.max(abs(wall_phi))),
                   capacitance_f=float(solution.capacitance))
        rows.append(row)
        np.savez_compressed(out/f'coefficients_panels_{n}.npz', starts=solution.sources.starts,
                            ends=solution.sources.ends, points=solution.sources.points, strengths=solution.strengths)
        print('cage', json.dumps(row), flush=True)
        (out/'refinement.json').write_text(json.dumps(rows, indent=1))
    solution.sources.quadrature *= 2
    high_quad = solution.evaluate(holdout)
    solution.sources.quadrature //= 2
    quadrature_error = float(np.max(abs(high_quad-solution.evaluate(holdout))))
    report = dict(inputs=inputs, reference=reference, refinement=rows,
                  common_off_fit_boundary_points=len(holdout),
                  independent_finite_line=line,
                  interior_vacuum_field=interior_field(solution, diameter/2),
                  image_quadrature_doubling_max_change_over_bias=quadrature_error,
                  limitations=['Vacuum only; no electron well or ion confinement is asserted.',
                               'Finite harmonic-source boundary fit: report off-fit errors, not just algebraic fit.',
                               'Spherical wall and rounded joints are assumptions; support and printed surface roughness absent.',
                               'No optimal bias, pressure, fuel, frequency or fusion gain has been calculated.'])
    (out/'report.json').write_text(json.dumps(report, indent=1, allow_nan=False))
    np.savez_compressed(out/'solution.npz', vertices=cage.vertices, edges=cage.edges,
                        source_starts=solution.sources.starts, source_ends=solution.sources.ends,
                        source_points=solution.sources.points, strengths=solution.strengths)
    plot(cage, solution, out)
    print(out/'report.json', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('out')
    parser.add_argument('--panels', nargs='+', type=int, default=[2, 4, 6])
    parser.add_argument('--diameter', type=float, default=.15)
    parser.add_argument('--wire', type=float, default=.002)
    parser.add_argument('--wall', type=float, default=.2)
    args = parser.parse_args()
    main(args.out, args.panels, args.diameter, args.wire, args.wall)
