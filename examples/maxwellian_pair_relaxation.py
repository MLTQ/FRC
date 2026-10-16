"""Complete collision relaxation compared with independent velocity integrals."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge as E, epsilon_0 as EPS0, electron_mass as ME, atomic_mass
from scipy.integrate import solve_ivp
from plasma_bench.maxwellian_pair_collisions import pair_moments
from .implicit_plasma_records import provenance, freeze


def velocity_integrals(state, density, masses, log, quadrature):
    """Integrate binary momentum/energy transfer over relative velocity space.

    Spherical polar integration around the mean drift; Gaussian conditional
    centre velocity performs the other three integrals. No erf or production
    collision function is used. Opposite polar directions are paired to retain
    small nonzero drift accurately. The tail is truncated 12 thermal speeds out.
    """
    ua, ub = state[:3], state[3:6]; ta, tb = E*state[6:]
    ma, mb = masses; total = ma+mb; reduced = ma*mb/total
    drift = ua-ub; speed = np.linalg.norm(drift)
    thermal = np.sqrt(ta/ma+tb/mb); d = speed/thermal
    xr, wr, xc, wc = quadrature
    y = .5*(d+12)*(xr+1); radial = .5*(d+12)*wr
    cosine = .5*(xc+1); angular = .5*wc
    yy = y[:, None]; cc = cosine[None, :]
    positive = np.exp(-.5*(yy-d)**2-yy*d*(1-cc))
    negative = np.exp(-.5*(yy-d)**2-yy*d*(1+cc))
    difference = positive*(-np.expm1(-2*yy*d*cc))
    total_gaussian = positive+negative
    weight = radial[:, None]*angular[None, :]/np.sqrt(2*np.pi)
    parallel = np.sum(weight*cc*difference)/thermal**2
    conditional = (ta-tb)/(total*thermal**2)
    centre_dot = conditional*(yy*total_gaussian-d*cc*difference)
    prefactor = 4*np.pi*np.prod(density)*(E**2/(4*np.pi*EPS0))**2*log/reduced
    force = -prefactor*parallel*drift/speed if speed else np.zeros(3)
    qa = -prefactor/thermal*np.sum(weight*(centre_dot-mb/total*d*cc*difference))
    qb = prefactor/thermal*np.sum(weight*(centre_dot+ma/total*d*cc*difference))
    return force, qa, qb


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh collision-relaxation directory required')
    out.mkdir(parents=True); identity = provenance(['examples.maxwellian_pair_relaxation']); freeze(out, identity)
    quadrature = (*np.polynomial.legendre.leggauss(96), *np.polynomial.legendre.leggauss(64))
    cases = [('equal_mass', (atomic_mass, atomic_mass), (2., 2.), (1e19, 1e19)),
             ('electron_argon', (ME, 39.948*atomic_mass), (10., 1.), (1e20, 1e20))]
    report = dict(source_sha256=identity['source_sha256'], terminal=False, cases={},
        scope='Homogeneous collision-only Maxwellian closure: no evolving fields, spatial transport or plasma-vacuum geometry. Independent velocity-space integral, not an experimental formation model.')
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for row, (name, masses, temperatures, density) in enumerate(cases):
        ma, mb = masses; na, nb = density; log = 10.
        thermal = np.sqrt(E*(temperatures[0]/ma+temperatures[1]/mb))
        relative = 3*thermal*np.array([.6, .8, 0.]); boost = np.array([200., -100., 300.])
        rhoa, rhob = na*ma, nb*mb
        initial = np.r_[boost+rhob/(rhoa+rhob)*relative,
                        boost-rhoa/(rhoa+rhob)*relative, temperatures]
        zero = pair_moments(na, nb, *temperatures, np.zeros(3), ma, mb, coulomb_log=log)
        timescale = 1/(float(zero['drag_coefficient'])*(1/rhoa+1/rhob))
        duration = 100*timescale

        def rhs(time, state, reference=False):
            if np.any(state[6:] <= 0) or not np.isfinite(state).all():
                raise ValueError('Nonpositive or nonfinite actual species state')
            if reference:
                force, qa, qb = velocity_integrals(state, density, masses, log, quadrature)
            else:
                a = pair_moments(na, nb, *state[6:], state[:3]-state[3:6], ma, mb, coulomb_log=log)
                force, qa, qb = a['force_on_a'], a['heat_to_a'], a['heat_to_b']
            return np.r_[force/rhoa, -force/rhob, qa/(1.5*na*E), qb/(1.5*nb*E)]

        scales = np.r_[np.full(6, thermal), temperatures]
        reference = solve_ivp(lambda t, scaled: rhs(t*timescale, scaled*scales, True)*timescale/scales,
            (0., duration/timescale), initial/scales, method='DOP853', rtol=2e-11, atol=2e-13, dense_output=True)
        if not reference.success:
            raise RuntimeError(reference.message)
        case = dict(density_m3=density, masses_kg=masses, initial_temperature_ev=temperatures,
            coulomb_log=log, initial_relative_speed_m_s=float(np.linalg.norm(relative)),
            duration_s=duration, initial_small_drift_relaxation_time_s=timescale, runs=[])
        energy0 = .5*rhoa*initial[:3]@initial[:3]+.5*rhob*initial[3:6]@initial[3:6]+1.5*E*(na*initial[6]+nb*initial[7])
        momentum0 = rhoa*initial[:3]+rhob*initial[3:6]
        equilibrium = (energy0-.5*(rhoa+rhob)*(boost@boost))/(1.5*E*(na+nb))
        case['conservation_equilibrium_temperature_ev'] = equilibrium
        for steps in (100, 200, 400):
            t = np.linspace(0., duration, steps+1); history = [initial.copy()]; state = initial.copy(); dt=duration/steps
            for k in range(steps):
                k1=rhs(t[k], state); k2=rhs(t[k]+dt/2, state+dt*k1/2)
                k3=rhs(t[k]+dt/2, state+dt*k2/2); k4=rhs(t[k+1], state+dt*k3)
                state=state+dt*(k1+2*k2+2*k3+k4)/6
                if np.any(state[6:] <= 0) or not np.isfinite(state).all():
                    raise ValueError('Rejected collision state; no thermal correction')
                history.append(state.copy())
            h = np.array(history); expected = reference.sol(t/timescale).T*scales
            energy = .5*rhoa*np.sum(h[:, :3]**2, axis=1)+.5*rhob*np.sum(h[:, 3:6]**2, axis=1)+1.5*E*(na*h[:, 6]+nb*h[:, 7])
            momentum = rhoa*h[:, :3]+rhob*h[:, 3:6]
            entropy = 1.5*(na*np.log(h[:, 6]/h[0, 6])+nb*np.log(h[:, 7]/h[0, 7]))
            verror = np.linalg.norm((h[:, :3]-h[:, 3:6])-(expected[:, :3]-expected[:, 3:6]), axis=1)
            terr = np.max(abs(h[:, 6:]-expected[:, 6:]), axis=0)
            case['runs'].append(dict(steps=steps, maximum_relative_drift_error=float(verror.max()/np.linalg.norm(relative)),
                maximum_temperature_errors_ev=terr.tolist(), maximum_relative_energy_drift=float(abs(energy/energy0-1).max()),
                maximum_momentum_error_over_initial_drift_momentum=float(np.linalg.norm(momentum-momentum0, axis=1).max()/(rhoa*np.linalg.norm(relative))),
                minimum_entropy_increment_over_total_particles=float(np.diff(entropy).min()/(na+nb)),
                final_temperature_ev=h[-1, 6:].tolist(), final_relative_speed_m_s=float(np.linalg.norm(h[-1,:3]-h[-1,3:6]))))
            np.savez_compressed(out/f'{name}_{steps}.npz', time_s=t, state=h, reference=expected,
                                physical_energy_j_m3=energy, entropy_change_over_kB_m3=entropy)
            axes[row, 0].plot(t/timescale, np.linalg.norm(h[:, :3]-h[:, 3:6], axis=1)/np.linalg.norm(relative), label=f'{steps} steps')
            axes[row, 2].plot(t/timescale, (energy-energy0)/energy0, label=f'{steps} steps')
            print(name, steps, case['runs'][-1], flush=True)
        axes[row, 0].plot(t/timescale, np.linalg.norm(expected[:,:3]-expected[:,3:6],axis=1)/np.linalg.norm(relative),'k--',label='Velocity-integral reference')
        axes[row, 1].plot(t/timescale, h[:,6],label='Species a');axes[row, 1].plot(t/timescale,h[:,7],label='Species b')
        axes[row, 1].plot(t/timescale,expected[:,6:],'k--',lw=.7)
        axes[row, 1].axhline(equilibrium,color='grey',ls=':',label='Eventual energy equilibrium')
        report['cases'][name]=case
        for ax,title,ylabel in zip(axes[row], ('Relative species drift','Species temperatures','Observed energy drift'), ('|Ua-Ub| / initial','Temperature [eV]','(E-E0) / E0')):
            ax.set(xlabel='Time / initial small-drift relaxation time',ylabel=ylabel,title=f'{name}: {title}');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Finite-drift Coulomb collision relaxation\nPaired drag and species heat; declared lnΛ=10; isotropic Maxwellian closure')
    fig.savefig(out/'relaxation.png',dpi=170);fig.savefig(out/'relaxation.pdf');plt.close(fig)
    if provenance(['examples.maxwellian_pair_relaxation']) != identity:
        raise RuntimeError('Collision source changed during run')
    report['terminal']=True
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination');args=parser.parse_args();run(args.destination)
