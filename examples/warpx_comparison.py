"""Run one isolated WarpX verification case; requires the optional GPU build."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import c, epsilon_0, mu_0, proton_mass
from .backend_wave_reference import ParallelWave, mode_amplitude, relative_error, vacuum_field


def run(args):
    # WarpX owns process-global geometry; every case is a separate process.
    from pywarpx import picmi, particle_containers
    case = ParallelWave(helicity=args.helicity, amplitude=args.amplitude)
    n = args.cells
    length = case.length if args.kind == 'hall' else 1.
    stop = 2/case.omega if args.kind == 'hall' else length/c
    limit = args.dt_omega/case.omega_i if args.kind == 'hall' else args.cfl*length/(n*c)
    count = int(np.ceil(stop/limit))
    dt = stop/count
    directory = Path(args.output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    grid = picmi.Cartesian1DGrid(number_of_cells=[n], lower_bound=[0.],
        upper_bound=[length], lower_boundary_conditions=['periodic'],
        upper_boundary_conditions=['periodic'], warpx_max_grid_size=n)
    if args.kind == 'hall':
        solver = picmi.HybridPICSolver(grid=grid, Te=0., n0=case.density,
            gamma=1., n_floor=args.floor_fraction*case.density,
            plasma_resistivity=0., plasma_hyper_resistivity=0., substeps=10)
    else:
        solver = picmi.ElectromagneticSolver(grid=grid, method='Yee')
    sim = picmi.Simulation(solver=solver, time_step_size=dt, max_steps=count,
        verbose=0, particle_shape=1, warpx_use_filter=False,
        warpx_amrex_the_arena_init_size=64*1024**2,
        warpx_grid_type='collocated' if args.kind == 'hall' else 'staggered',
        warpx_used_inputs_file=str(directory/'warpx_used_inputs'),
        warpx_current_deposition_algo='direct' if args.kind == 'hall' else 'esirkepov')
    if args.kind == 'hall':
        bx = f'{case.amplitude:.17g}*cos({case.k:.17g}*z)'
        by = f'{-args.helicity*case.amplitude:.17g}*sin({case.k:.17g}*z)'
        ions = picmi.Species(name='ions', charge=picmi.constants.q_e,
            mass=proton_mass, initial_distribution=picmi.AnalyticDistribution(
                density_expression=f'{case.density:.17g}',
                momentum_expressions=[f'{case.velocity_ratio:.17g}*({bx})',
                                      f'{case.velocity_ratio:.17g}*({by})', '0']))
        sim.add_species(ions, layout=picmi.GriddedLayout(
            grid=grid, n_macroparticle_per_cell=[args.ppc]))
        initial = picmi.AnalyticInitialField(Bx_expression=bx, By_expression=by,
            Bz_expression=case.bias, Ex_expression=f'{case.omega/case.k:.17g}*({by})',
            Ey_expression=f'{-case.omega/case.k:.17g}*({bx})', Ez_expression='0')
    else:
        by = f'{args.amplitude:.17g}*cos(2*pi*z/{length:.17g})'
        initial = picmi.AnalyticInitialField(Bx_expression='0', By_expression=by,
            Bz_expression='0', Ex_expression=f'{c:.17g}*({by})',
            Ey_expression='0', Ez_expression='0')
    sim.add_applied_field(initial)
    sim.add_diagnostic(picmi.ReducedDiagnostic(name='field_energy',
        diag_type='FieldEnergy', period=1, path=str(directory)+'/', extension='txt'))
    if args.kind == 'hall':
        sim.add_diagnostic(picmi.ReducedDiagnostic(name='particle_energy',
            diag_type='ParticleEnergy', period=1, path=str(directory)+'/', extension='txt'))
    begin = time.perf_counter()
    sim.initialize_inputs()
    sim.initialize_warpx()
    initialization_s = time.perf_counter()-begin

    def read(component, axis):
        mf = sim.fields.get(component+'field_fp', dir=axis, level=0)
        # Public global indexing gathers device data to host and drops ghosts.
        data = np.asarray(mf[:]).squeeze().copy()
        if data.size not in (n, n+1):
            raise RuntimeError(f'Unexpected field shape: {data.shape}')
        # Nodal components include the periodic duplicate; cell centers do not.
        offset = 0. if data.size == n+1 else .5
        return data[:n], (np.arange(n)+offset)*length/n

    axes = [('B', 'x'), ('B', 'y'), ('B', 'z'), ('E', 'x'), ('E', 'y'), ('E', 'z')]
    samples = {}
    def particles(state):
        if args.kind != 'hall':
            return
        container = particle_containers.ParticleContainerWrapper('ions')
        for axis in ('ux', 'uy', 'uz', 'weight'):
            arrays = getattr(container, 'get_particle_'+axis)(copy_to_host=True)
            samples['particle_'+axis+'_'+state] = np.concatenate(arrays)

    particles('initial')
    for component, axis in axes:
        values, coordinates = read(component, axis)
        samples[component+axis+'_initial'] = values
        samples[component+axis+'_z'] = coordinates
    begin = time.perf_counter()
    sim.step(count)
    actual_time = float(sim.extension.warpx.gett_new(lev=0))
    actual_steps = int(sim.extension.warpx.getistep(lev=0))
    if actual_steps != count or not np.isclose(actual_time, stop, rtol=1e-12, atol=0.):
        raise RuntimeError(f'Incomplete evolution: {actual_steps} steps, time {actual_time}')
    for component, axis in axes:
        values, _ = read(component, axis)
        samples[component+axis+'_final'] = values
    particles('final')
    elapsed = time.perf_counter()-begin  # Includes endpoint readback/synchronization.
    result = {'completed': True, 'engine': 'WarpX', 'kind': args.kind,
        'arguments': vars(args),
        'cells': n, 'steps': count, 'dt_s': dt, 'stop_s': stop,
        'actual_time_s': actual_time, 'actual_steps': actual_steps,
        'initialization_s': initialization_s, 'evolution_and_readback_s': elapsed,
        'parameters': asdict(case) if args.kind == 'hall' else {'amplitude_T': args.amplitude},
        'numerics': {'particle_shape': 1, 'particles_per_cell': args.ppc if args.kind == 'hall' else 0,
                     'filter': False, 'resistivity': 0., 'hyper_resistivity': 0.,
                     'density_floor_m-3': args.floor_fraction*case.density if args.kind == 'hall' else None,
                     'hybrid_field_substeps': 10 if args.kind == 'hall' else None}}
    if args.kind == 'hall':
        z = samples['Bx_z']
        actual = np.array([samples['Bx_final'], samples['By_final']])
        exact = case.magnetic(z, stop)
        amplitude = mode_amplitude(samples['Bx_final'], case.k*z)/case.amplitude
        result.update(relative_magnetic_error=relative_error(actual, exact),
            phase_error_rad=float(np.angle(amplitude*np.exp(1j*case.omega*stop))),
            amplitude_error=float(abs(amplitude)-1), omega_s_1=case.omega,
            initial_magnetic_error=relative_error(np.array([samples['Bx_initial'],
                samples['By_initial']]), case.magnetic(z, 0.)))
        samples['Bperp_exact_final'] = exact
        energies = []
        for state in ('initial', 'final'):
            u2 = sum(samples['particle_u'+a+'_'+state]**2 for a in 'xyz')
            kinetic = np.sum(proton_mass*samples['particle_weight_'+state]
                             *u2/(np.sqrt(1+u2/c**2)+1))
            # Remove unchanged bias energy so conservation is not diluted by B0.
            magnetic = length*sum(np.mean(samples['B'+a+'_'+state]**2)
                                  for a in 'xy')/(2*mu_0)
            energies.append(float(kinetic+magnetic))
        result['perturbation_energy_J_m-2'] = energies
        result['relative_endpoint_energy_change'] = (energies[1]-energies[0])/energies[0]
    else:
        for component in ('By', 'Ex'):
            exact = vacuum_field(samples[component+'_z'], stop, length,
                                 args.amplitude, electric=component == 'Ex')
            result[component+'_relative_error'] = relative_error(samples[component+'_final'], exact)
            samples[component+'_exact_final'] = exact
        energies = []
        for state in ('initial', 'final'):
            energies.append(sum(np.mean(samples['B'+a+'_'+state]**2)/(2*mu_0)
                + epsilon_0*np.mean(samples['E'+a+'_'+state]**2)/2 for a in 'xyz')*length)
        result['relative_endpoint_energy_change'] = (energies[1]-energies[0])/energies[0]
    names = [Path(__file__), Path(__file__).with_name('backend_wave_reference.py')]
    result['source_files'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in names}
    archive = directory/'wave_samples.npz'
    np.savez_compressed(archive, **samples)
    result['archive_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
    (directory/'measurement.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
    sim.finalize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=['hall', 'vacuum'], required=True)
    parser.add_argument('--cells', type=int, default=32)
    parser.add_argument('--helicity', type=int, choices=[-1, 1], default=1)
    parser.add_argument('--ppc', type=int, default=16)
    parser.add_argument('--dt-omega', type=float, default=.005)
    parser.add_argument('--floor-fraction', type=float, default=1e-8)
    parser.add_argument('--cfl', type=float, default=.5)
    parser.add_argument('--amplitude', type=float, default=1e-7)
    parser.add_argument('--output', required=True)
    run(parser.parse_args())
