"""Independent annular-vacuum and absorbing-wall experiments (not physical FRC formation)."""
import json
from pathlib import Path
import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0 as MU0, m_p as MP
from plasma_bench.hybrid3d import curl_edges, curl_faces, divergence_faces, HybridParams
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, GapProjector, DeviceHybrid, GAP


def ring_projection(dx):
    """Exact Biot–Savart A of a 10 kA circular filament inside the gas annulus.

    Mapping its curl makes compatible boundary-normal data; the vacuum result
    is compared with the independent elliptic-integral B, not with a second
    copy of the projector. No particle/fluid approximation enters the source.
    """
    dev = Device(.12, .12, .15, .03, .04)
    box = device_box(dev, dx, dx)
    types = cell_types(box, dev)
    projector = GapProjector(types, box.d)
    xs = [box.lo[a]+box.d[a]*np.arange(box.shape[a]) for a in range(3)]
    xh = [x+.5*box.d[a] for a, x in enumerate(xs)]
    radius, current = .08, 1e4

    def analytic(x, y, z):
        r = np.hypot(x, y)
        alpha, beta = (radius-r)**2+z*z, (radius+r)**2+z*z
        m = np.clip(4*radius*r/beta, 1e-20, 1-1e-12)
        k, e = ellipk(m), ellipe(m)
        fac = MU0*current/(2*np.pi*np.sqrt(beta))
        bz = fac*(k+(radius**2-r*r-z*z)/np.maximum(alpha, 1e-20)*e)
        br = fac*z/np.maximum(r, 1e-30)*(-k+(radius**2+r*r+z*z)/np.maximum(alpha, 1e-20)*e)
        br = np.where(r < 1e-10, 0., br)
        ap = MU0*current/np.pi*np.sqrt(radius/np.maximum(r*m, 1e-30))*((1-m/2)*k-e)
        ap = np.where(r < 1e-10, 0., ap)
        return r, br, bz, ap
    edges = []
    for a in range(2):
        x, y, z = np.meshgrid(*[xh[q] if q == a else xs[q] for q in range(3)], indexing='ij')
        r, _, _, ap = analytic(x, y, z)
        edges.append(ap*(-y if a == 0 else x)/np.maximum(r, 1e-30))
    b = list(curl_edges(*edges, np.zeros(box.shape), box.d))
    exact, selection = [], []
    for a in range(3):
        x, y, z = np.meshgrid(*[xs[q] if q == a else xh[q] for q in range(3)], indexing='ij')
        r, br, bz, _ = analytic(x, y, z)
        exact.append(bz if a == 2 else br*(x if a == 0 else y)/np.maximum(r, 1e-30))
        selection.append((types == GAP) & (np.roll(types, 1, a) == GAP))
    plain = projector([v.copy() for v in b])
    projector.meridional_harmonic(box, dev.inner_radius)
    projected = projector([v.copy() for v in b])
    norm = sum(np.sum(exact[a][selection[a]]**2) for a in range(3))
    error = lambda field: float(np.sqrt(sum(np.sum((field[a][selection[a]]-exact[a][selection[a]])**2)
                                               for a in range(3))/norm))
    h, _ = projector.harmonics[0]
    hj = curl_faces(*h, box.d)
    interior_curl = []
    for a in range(3):
        axes = [q for q in range(3) if q != a]
        g = types == GAP
        allgap = g & np.roll(g, 1, axes[0]) & np.roll(g, 1, axes[1]) & np.roll(np.roll(g, 1, axes[0]), 1, axes[1])
        interior_curl.append(float(np.abs(hj[a][allgap]).max()))
    return dict(dx_m=dx, relative_B_l2=error(projected), error_without_harmonic=error(plain),
                max_gap_div_B=float(np.abs(divergence_faces(*projected, box.d)[types == GAP]).max()),
                harmonic_curl_max=max(interior_curl))


def ballistic_wall(dx):
    """Cold ions fly inward with E=B=0; exact circular-wall impact is (r0-a)/v."""
    dev = Device(.12, .12, .15, .03, .04)
    box = device_box(dev, dx, dx)
    angle = (np.arange(128)+.37)*np.pi*(3-np.sqrt(5))
    direction = np.column_stack((np.cos(angle), np.sin(angle), np.zeros(len(angle))))
    x, v = .08*direction, -1e5*direction
    model = DeviceHybrid(box, dev, HybridParams(ion_mass=MP, n_floor=1.),
                         [np.zeros(box.shape) for _ in range(3)], x, v, 1., .25*dx/1e5, n_ref=1.)
    loss_times = []
    while model.x.shape[0] and model.time < 1e-6:
        previous = len(model.x)
        model.step()
        loss_times.extend([model.time]*(previous-len(model.x)))
    exact = (.08-.04)/1e5
    return dict(dx_m=dx, expected_impact_s=exact, impact_rms_error_s=float(np.sqrt(np.mean((np.array(loss_times)-exact)**2))),
                count_balance=model.particles()+model.lost['count'], lost=model.lost,
                expected_absorbed_energy_j=.5*MP*128*1e10)


def driven_vacuum():
    """Empty chamber: an arbitrarily tilted uniform B ramp is an exact vacuum solution."""
    dev = Device(0., .12, .15, .03)
    box = device_box(dev, .02, .02)
    initial = np.array([.08, -.03, .12])
    duration, slope = 4e-7, .2/4e-7

    class Ramp:
        def __call__(self, time_s):
            return np.array([np.full(box.shape, b*(1+slope*time_s)) for b in initial])

        def edge_derivative(self, time_s):
            fields = []
            for a in range(3):
                xs = [box.lo[q]+box.d[q]*(np.arange(box.shape[q])+(.5 if a == q else 0.)) for q in range(3)]
                x, y, z = np.meshgrid(*xs, indexing='ij')
                cross = [initial[1]*z-initial[2]*y, initial[2]*x-initial[0]*z, initial[0]*y-initial[1]*x]
                fields.append(.5*slope*cross[a])
            return fields
    source = Ramp()
    model = DeviceHybrid(box, dev, HybridParams(ion_mass=MP, n_floor=1.), source(0.),
                         np.empty((0, 3)), np.empty((0, 3)), 1., 2e-8, n_ref=1., source=source)
    maximum = 0.
    gap = model.types == GAP
    while model.time < duration-1e-15:
        model.step()
        exact = source(model.time)
        maximum = max(maximum, max(float(np.abs(model.b[a][gap]-exact[a][gap]).max()) for a in range(3)))
    return dict(steps=model.steps, max_B_error_T=maximum, max_relative_div_B=model.max_divergence(),
                final_B_T=(initial*(1+slope*model.time)).tolist(), particles=model.particles())


def main(out='outputs/hybrid3d_annular_check'):
    out = Path(out); out.mkdir(exist_ok=True, parents=True)
    result = dict(interpretation='Verification fixtures, not FRC experiments',
                  ring=[ring_projection(dx) for dx in (.02, .01, .005)],
                  ballistic=[ballistic_wall(dx) for dx in (.02, .01, .005)],
                  driven_vacuum=driven_vacuum())
    (out/'comparison.json').write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1), flush=True)
    return result


if __name__ == '__main__':
    main()
