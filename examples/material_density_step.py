"""Complete strong-density-jump Euler evolution with an exact central reference."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as EV
from scipy.optimize import brentq
from plasma_bench.precise_geometry import coordinates, quadrature
from plasma_bench.precise_material import PreciseMaterialSubcells, PreciseMaterialSet
from plasma_bench.material_species_set import MaterialSetState
from .precise_material import save_states
from .implicit_plasma_records import provenance, freeze


def exact_jump(position, time, ratio=1e-8, gamma=5/3):
    """Zero initial flow, equal temperature; normalized left rho=p=1."""
    if time <= 0 or not 0 < ratio < 1:
        raise ValueError('Positive time and a pressure/density ratio strictly between zero and one required')
    def wave(p, rho, pk):
        sound = np.sqrt(gamma*pk/rho)
        if p > pk:
            return (p-pk)*np.sqrt(2/((gamma+1)*rho)/(p+(gamma-1)/(gamma+1)*pk))
        return 2*sound/(gamma-1)*np.expm1((gamma-1)/(2*gamma)*np.log(p/pk))
    logp = brentq(lambda y: wave(np.exp(y), 1., 1.)+wave(np.exp(y), ratio, ratio),
                  np.log(ratio), 0., xtol=2e-14, rtol=2e-14)
    pstar = np.exp(logp)
    speed = .5*(wave(pstar, ratio, ratio)-wave(pstar, 1., 1.))
    sound = np.sqrt(gamma)
    star_sound = sound*pstar**((gamma-1)/(2*gamma))
    shock = sound*np.sqrt((gamma+1)/(2*gamma)*pstar/ratio+(gamma-1)/(2*gamma))
    xi = np.asarray(position, float)/time
    rho, v, p = np.full_like(xi, ratio), np.zeros_like(xi), np.full_like(xi, ratio)
    left = xi < -sound
    fan = (xi >= -sound) & (xi < speed-star_sound)
    left_star = (xi >= speed-star_sound) & (xi < speed)
    right_star = (xi >= speed) & (xi < shock)
    rho[left], p[left] = 1., 1.
    cs = 2/(gamma+1)*(sound-(gamma-1)*xi[fan]/2)
    rho[fan], p[fan] = (cs/sound)**(2/(gamma-1)), (cs/sound)**(2*gamma/(gamma-1))
    v[fan] = 2/(gamma+1)*(sound+xi[fan])
    rho[left_star], p[left_star], v[left_star] = pstar**(1/gamma), pstar, speed
    beta = (gamma-1)/(gamma+1)
    rho_right_star = ratio*(pstar/ratio+beta)/(beta*pstar/ratio+1)
    rho[right_star], p[right_star], v[right_star] = rho_right_star, pstar, speed
    details = dict(star_pressure=pstar, contact_speed=speed,
        wave_speeds=[-sound, speed-star_sound, speed, shock],
        maximum_temperature_ratio=max(1., pstar/rho_right_star),
        pressure_root_residual=wave(pstar, 1., 1.)+wave(pstar, ratio, ratio))
    return (rho, v, p, p/rho), details


def run(destination):
    out = Path(destination)
    if out.exists():
        raise ValueError('Fresh density-step output required')
    out.mkdir(parents=True)
    identity = provenance(['examples.material_density_step'])
    freeze(out, identity)
    began = time.perf_counter()
    mass, temperature, density, a, tau = 39.948*atomic_mass, .2, 1e14, .01, .08
    speed0 = np.sqrt(EV*temperature/mass)
    duration = tau*a/speed0
    _, details = exact_jump([0.], tau)
    if max(abs(np.array(details['wave_speeds'])))*tau >= .5 or 1-np.sqrt(5/3)*tau <= .5:
        raise ValueError('Wave interactions invalidate the stated central comparison')
    report = dict(terminal=False, source_sha256=identity['source_sha256'], requested_time_s=duration,
        initial_temperature_ev=temperature, initial_left_density_m3=density, right_to_left_density_ratio=1e-8,
        half_width_m=a, reference=details, transport_fraction=.5, cases=[],
        scope='Complete neutral Euler pressure/density-step fixture with physical argon mass and exact external vacuum. '
              'Independent central rarefaction/contact/shock before exterior fans arrive; no reacting/FRC validation.')
    def save():
        report['elapsed_seconds'] = time.perf_counter()-began
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    xplot = np.linspace(-.5, .5, 3001)
    reference, _ = exact_jump(xplot, tau)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    for ax, y in zip(axes, (reference[1]*speed0, reference[3]*temperature, np.log10(reference[0]*density))):
        ax.plot(xplot*a*1e3, y, 'k--', label='Exact central Euler solution')
    for cells in (64, 128, 256):
        faces = coordinates(np.linspace(-a, a, cells+1))
        widths = np.array([float(c-b) for b, c in zip(faces[:-1], faces[1:])])
        number = density*np.where(np.asarray(faces[:-1], float) < 0, 1., 1e-8)*widths
        component = PreciseMaterialSubcells(number, mass, 0., (-.2, .2), faces, cells,
            velocity_faces=faces, numerical_viscosity=(.5, 2.), viscosity_limiter='minmod')
        model = PreciseMaterialSet([component], [0])
        state = MaterialSetState((component.state(faces[None, :], temperature),))
        initial = model.totals(state)
        times, states, energies, momenta = [0.], [state], [initial['energy_j_m2']], [initial['momentum_kg_per_m_s']]
        failure = None
        while times[-1] < duration:
            dt = min(.5*model.timestep(state), duration-times[-1])
            try:
                if not np.isfinite(dt) or dt <= 0 or times[-1]+dt == times[-1]:
                    raise ValueError('Physical transport clock cannot advance')
                if len(times) > 100000:
                    raise ValueError('Declared100000-step limit; full duration not complete')
                state = model.step(state, dt)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
                failure = str(error)
                break
            times.append(min(times[-1]+dt, duration))
            states.append(state)
            totals = model.totals(state)
            energies.append(totals['energy_j_m2'])
            momenta.append(totals['momentum_kg_per_m_s'])
        breaks = coordinates(np.array(details['wave_speeds'])*times[-1]*speed0)
        numerical_faces = state.blocks[0].faces[0]
        edges = np.unique(np.r_[coordinates([-a/2, a/2]), breaks,
            numerical_faces[(numerical_faces > -a/2) & (numerical_faces < a/2)]])
        points, weights = quadrature(edges, 16)
        x, w = points.ravel(), weights.ravel()/a
        sampled = model.sample(state, x)
        actual = (sampled['number'][0]/density, sampled['velocity'][0, :, 2]/speed0,
                  sampled['number'][0]*sampled['temperature_ev'][0]/(density*temperature),
                  sampled['temperature_ev'][0]/temperature)
        if times[-1] > 0:
            exact, _ = exact_jump(np.asarray(x, float)/a, times[-1]*speed0/a)
        else:
            initial_density = np.where(np.asarray(x, float) < 0, 1., 1e-8)
            exact = (initial_density, np.zeros(len(x)), initial_density, np.ones(len(x)))
        errors = {key:float(np.sum(w*abs(y-ref))) for key, y, ref in zip(
            ('density', 'velocity', 'pressure', 'temperature'), actual, exact)}
        energy = np.asarray(energies)
        row = dict(cells=cells, success=bool(failure is None and times[-1] == duration), failure=failure,
            completed_time_s=times[-1], steps=len(times)-1, central_normalized_l1=errors,
            maximum_central_temperature_ev=float(sampled['temperature_ev'][0].max()),
            exact_maximum_central_temperature_ev=temperature*details['maximum_temperature_ratio'],
            maximum_relative_energy_error=float(abs(energy-energy[0]).max()/energy[0]),
            maximum_momentum_change_kg_per_m_s=float(np.linalg.norm(np.asarray(momenta)-momenta[0], axis=1).max()))
        save_states(out/f'n{cells}_states.npz', [model]*len(states), states)
        np.savez_compressed(out/f'n{cells}_budget.npz', time_s=times, dt_s=np.diff(times),
            energy_j_m2=energies, momentum_kg_per_m_s=momenta)
        report['cases'].append(row)
        save()
        print(row, flush=True)
        sampled = model.sample(state, coordinates(xplot*a))
        for ax, y in zip(axes, (sampled['velocity'][0, :, 2], sampled['temperature_ev'][0], np.log10(sampled['number'][0]))):
            ax.plot(xplot*a*1e3, y, label=f'{cells} velocity cells', alpha=.8)
    for ax, label in zip(axes, ('Velocity [m/s]', 'Temperature [eV]', 'log₁₀(number density / m⁻³)')):
        ax.set(xlabel='Axial position [mm]', ylabel=label)
        ax.grid(alpha=.2)
    fig.legend(*axes[0].get_legend_handles_labels(), loc='outside lower center', ncol=4, fontsize=8)
    fig.suptitle('Equal-temperature gas expanding across a 10⁸ density contrast')
    fig.savefig(out/'density_step.png', dpi=160)
    plt.close(fig)
    if identity != provenance(['examples.material_density_step']):
        raise RuntimeError('Source changed during density-step evolution')
    report['terminal'] = True
    report['success'] = all(row['success'] for row in report['cases'])
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('destination')
    args = parser.parse_args()
    result = run(args.destination)
    raise SystemExit(0 if result['success'] else 1)
