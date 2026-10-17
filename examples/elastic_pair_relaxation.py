"""Complete neutral elastic relaxation with independent binary velocity integrals."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as EV, electron_mass as ME, atomic_mass
from scipy.integrate import solve_ivp
from plasma_bench.elastic_maxwellian import TabulatedElasticPair
from plasma_bench.lxcat_cross_sections import CrossSection, read_cross_sections
from .implicit_plasma_records import provenance, freeze


def velocity_integrals(state, density, masses, section, quadrature):
    """Independent radial/polar binary-event integral; no analytic angular kernel."""
    ma, mb = masses; total = ma+mb; reduced = ma*mb/total
    ta, tb = EV*state[6:]; u = state[:3]-state[3:6]
    speed = np.linalg.norm(u); variance = ta/ma+tb/mb; thermal = np.sqrt(variance)
    xr, wr, xc, wc = quadrature
    edges = np.sqrt(2*EV*section.energy_ev/ma)
    added = speed+thermal*np.arange(-12., 13.)
    edges = np.unique(np.r_[edges, added[(added > edges[0]) & (added < edges[-1])]])
    width = np.diff(edges)
    g = (edges[:-1, None]+.5*width[:, None]*(1+xr)).ravel()
    radial = (.5*width[:, None]*wr/thermal).ravel()
    r, d = g[:, None]/thermal, speed/thermal
    c, angular = .5*(xc+1)[None, :], .5*wc[None, :]
    positive = np.exp(-.5*(r-d)**2-r*d*(1-c))
    negative = np.exp(-.5*(r-d)**2-r*d*(1+c))
    difference = positive*(-np.expm1(-2*r*d*c))
    combined = positive+negative
    sigma = np.interp(.5*ma*g*g/EV, section.energy_ev, section.sigma_m2)
    weight = radial[:, None]*angular*r*r/np.sqrt(2*np.pi)*sigma[:, None]
    parallel = np.sum(weight*g[:, None]**2*c*difference)
    centred = np.sum(weight*g[:, None]**2*(g[:, None]*combined-speed*c*difference))
    factor = reduced*np.prod(density)
    force = -factor*parallel*u/speed if speed else np.zeros(3)
    exchange = factor*(tb-ta)/(total*variance)*centred
    power = factor*speed*parallel
    return force, mb/total*power+exchange, ma/total*power-exchange


def run(destination, cross_sections):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh elastic-relaxation output directory required')
    data = read_cross_sections(cross_sections)
    elastic = [s for s in data['sections'] if s.kind == 'ELASTIC' and s.target == 'Ar']
    if len(elastic) != 1:
        raise ValueError('Exactly one explicit Ar elastic momentum-transfer table required')
    out.mkdir(parents=True)
    identity = provenance(['examples.elastic_pair_relaxation']); freeze(out, identity)
    synthetic = CrossSection('ELASTIC', 'synthetic constant', np.array([]), (),
                             np.array([0., 1e5]), np.array([2e-19, 2e-19]))
    cases = [('constant_equal_mass', (39.948*atomic_mass,)*2, (2., 1.), (1e19, 1e19), synthetic),
             ('electron_argon', (ME, 39.948*atomic_mass), (2., .02585), (1e18, 1e20), elastic[0])]
    report = dict(terminal=False, source_sha256=identity['source_sha256'],
        private_input_sha256=data['input_sha256'], provider='LXCat BSR-500, O. Zatsarinny and K. Bartschat',
        provider_url='https://www.lxcat.net/BSR', table_energy_bounds_ev=elastic[0].energy_ev[[0,-1]].tolist(),
        scope='Homogeneous fixed-population elastic source, no fields or reactions. Isotropic Maxwellian closure; finite supplied table domain, no extrapolation. Private tables and derived rate curves are not bundled.', cases={})
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout='constrained')
    for row, (name, masses, temperatures, density, section) in enumerate(cases):
        ma, mb = masses; na, nb = density; rhoa, rhob = na*ma, nb*mb
        model = TabulatedElasticPair(section, ma, mb)
        thermal = np.sqrt(EV*(temperatures[0]/ma+temperatures[1]/mb))
        relative = 3*thermal*np.array([.6, .8, 0.]); boost = np.array([200., -100., 300.])
        initial = np.r_[boost+rhob/(rhoa+rhob)*relative, boost-rhoa/(rhoa+rhob)*relative, temperatures]
        k0 = model.coefficients(*temperatures, 0.)['drag']
        timescale = 1/(k0*na*nb*(1/rhoa+1/rhob)); duration = 8*timescale
        case = dict(masses_kg=masses, density_m3=density, initial_state=initial.tolist(),
                    duration_s=duration, initial_small_drift_relaxation_time_s=timescale, runs=[])
        if name == 'constant_equal_mass':
            exact = 8/3*np.sqrt(2/np.pi)*(ma*mb/(ma+mb))*2e-19*thermal
            case['analytic_zero_drift_drag_relative_error'] = abs(k0/exact-1)

        def rhs(state, quadrature=None):
            if not np.isfinite(state).all() or np.any(state[6:] <= 0):
                raise ValueError('Nonpositive actual temperature; no correction')
            if quadrature is None:
                value = model.moments(na, nb, *state[6:], state[:3]-state[3:6])
                force, qa, qb = value['force_on_a'], value['heat_to_a'], value['heat_to_b']
            else:
                force, qa, qb = velocity_integrals(state, density, masses, section, quadrature)
            return np.r_[force/rhoa, -force/rhob, qa/(1.5*na*EV), qb/(1.5*nb*EV)]

        scales = np.r_[np.full(6, thermal), temperatures]
        references = []
        for radial, angular in [(8, 32), (16, 64)]:
            q = (*np.polynomial.legendre.leggauss(radial), *np.polynomial.legendre.leggauss(angular))
            ref = solve_ivp(lambda t, y: rhs(y*scales, q)*timescale/scales,
                            (0., 8.), initial/scales, method='DOP853', rtol=2e-11, atol=2e-13, dense_output=True)
            if not ref.success:
                raise RuntimeError(ref.message)
            references.append(ref)
            print(name, 'independent reference complete', radial, angular, flush=True)
        check = np.linspace(0., 8., 401)
        reference_delta = (references[0].sol(check)-references[1].sol(check)).T*scales
        case['reference_quadrature_maximum_velocity_difference_m_s'] = float(abs(reference_delta[:, :6]).max())
        case['reference_quadrature_maximum_temperature_difference_ev'] = float(abs(reference_delta[:, 6:]).max())
        energy0 = .5*rhoa*initial[:3]@initial[:3]+.5*rhob*initial[3:6]@initial[3:6]+1.5*EV*(na*initial[6]+nb*initial[7])
        momentum0 = rhoa*initial[:3]+rhob*initial[3:6]
        case['conservation_equilibrium_temperature_ev'] = (energy0-.5*(rhoa+rhob)*(boost@boost))/(1.5*EV*(na+nb))
        for steps in (100, 200, 400):
            t = np.linspace(0., duration, steps+1); dt = duration/steps
            state = initial.copy(); history = [state.copy()]
            for _ in range(steps):
                a = rhs(state); b = rhs(state+dt*a/2); c = rhs(state+dt*b/2); d = rhs(state+dt*c)
                state = state+dt*(a+2*b+2*c+d)/6
                if not np.isfinite(state).all() or np.any(state[6:] <= 0):
                    raise ValueError('Rejected complete-step state; no thermal correction')
                history.append(state.copy())
            h = np.array(history); expected = references[-1].sol(t/timescale).T*scales
            energy = .5*rhoa*np.sum(h[:, :3]**2, axis=1)+.5*rhob*np.sum(h[:, 3:6]**2, axis=1)+1.5*EV*(na*h[:, 6]+nb*h[:, 7])
            momentum = rhoa*h[:, :3]+rhob*h[:, 3:6]
            error = (h[:, :3]-h[:, 3:6])-(expected[:, :3]-expected[:, 3:6])
            outside = [model.coefficients(*s[6:], np.linalg.norm(s[:3]-s[3:6]))['probability_outside_table'] for s in h]
            measurement = dict(steps=steps, maximum_drift_error_over_initial=float(np.linalg.norm(error, axis=1).max()/np.linalg.norm(relative)),
                maximum_temperature_errors_ev=np.max(abs(h[:, 6:]-expected[:, 6:]), axis=0).tolist(),
                maximum_relative_energy_drift=float(abs(energy/energy0-1).max()),
                maximum_scaled_momentum_error=float(np.linalg.norm(momentum-momentum0, axis=1).max()/(rhoa*np.linalg.norm(relative))),
                maximum_relative_probability_outside_table=float(max(outside)), final_temperature_ev=h[-1, 6:].tolist(),
                final_relative_speed_m_s=float(np.linalg.norm(h[-1, :3]-h[-1, 3:6])))
            case['runs'].append(measurement)
            np.savez_compressed(out/f'{name}_{steps}.npz', time_s=t, state=h, reference=expected, physical_energy_j_m3=energy)
            print(name, measurement, flush=True)
            axes[row, 0].plot(t/timescale, np.linalg.norm(h[:, :3]-h[:, 3:6], axis=1)/np.linalg.norm(relative), label=f'{steps} steps')
            axes[row, 2].plot(t/timescale, energy/energy0-1, label=f'{steps} steps')
        axes[row, 0].plot(t/timescale, np.linalg.norm(expected[:, :3]-expected[:, 3:6], axis=1)/np.linalg.norm(relative), 'k--', label='Velocity integral')
        axes[row, 1].plot(t/timescale, h[:, 6], label='Projectile'); axes[row, 1].plot(t/timescale, h[:, 7], label='Target')
        axes[row, 1].plot(t/timescale, expected[:, 6:], 'k--', lw=.7)
        for ax, title, label in zip(axes[row], ['Relative drift', 'Temperature', 'Energy drift'], ['Drift / initial', 'Temperature [eV]', '(E-E0) / E0']):
            ax.set(title=f'{name}: {title}', xlabel='Time / initial relaxation time', ylabel=label); ax.grid(alpha=.2); ax.legend(fontsize=7)
        report['cases'][name] = case
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    fig.suptitle('Elastic Maxwellian source relaxation: paired drag, recoil and heat\nComplete fixed-population components; no electromagnetic drive or reactions')
    fig.savefig(out/'relaxation.png', dpi=160); fig.savefig(out/'relaxation.pdf'); plt.close(fig)
    if provenance(['examples.elastic_pair_relaxation']) != identity:
        raise RuntimeError('Source changed during calculation')
    report['terminal'] = True
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('--cross-sections', required=True)
    args = parser.parse_args(); run(args.destination, args.cross_sections)
