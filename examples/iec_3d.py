"""Local, source-frozen electrostatic PIC experiment; no prescribed ion well."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.constants import elementary_charge as qe, epsilon_0, electron_mass

from plasma_bench.iec_field import Field
from plasma_bench.iec_mesh import build_mesh
from plasma_bench.iec_particles import WallEmitter, load_plasma
from plasma_bench.hops_drive import CageDrive
from plasma_bench.iec_diagnostics import ion_motion, ion_sensor
from plasma_bench.iec_checkpoint import write_checkpoint, load_checkpoint, digest, physics_hashes as source_hashes


DEFAULTS = dict(bias_v=100., modulation=0., frequency_mhz=1., rise_ns=5.,
                duration_ns=50., dt_ps=20., particles=600, density_m3=1e11,
                seed_radius_m=.05, electron_ev=1., ion_ev=.05, ion_amu=2.,
                source_ua=100., source_ev=2., mesh_mm=18., wire_mesh_mm=1., seed=41,
                drive_mode='legacy', response_ns=20., sensor_ns=10., warmup_ns=200.,
                min_dwell_ns=20., release_fraction=.8, return_fraction=1.,
                min_control_ions=30, pulse_period_ns=1000., pulse_duty=.5,
                saved_intervals=100)


def atomic_json(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, allow_nan=False, separators=(',', ':')))
    temporary.replace(path)


def voltage(config, t):
    ramp = 1. if config['rise_ns'] == 0 else 1-np.exp(-t/(config['rise_ns']*1e-9))
    return config['bias_v']*ramp*(1+config['modulation']*np.sin(2*np.pi*config['frequency_mhz']*1e6*t))


def run(output, config, mesh_path=None, *, checkpoint=None, resume=False, stop_step=None):
    output = Path(output).resolve()
    if resume and checkpoint is None:
        raise ValueError('Resume requires a complete checkpoint')
    if not resume:
        output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    if resume:
        state=json.loads((output/'run.json').read_text())
        if state['config'] != config:
            raise ValueError('Cannot alter an existing experiment on resume')
        restart=output/'restarts'/f'{time.time_ns()}_{os.getpid()}'
        restart.mkdir(parents=True)
        shutil.copy2(output/'run.json',restart/'previous_run.json')
        state.pop('error',None)
        state.update(status='restarting',restart_checkpoint=str(Path(checkpoint).resolve()))
    else:
        state = dict(id=output.name, experiment='iec', label=config.get('run_label', 'Positive cage · 3D kinetic'),
                     status='meshing', config=config, frames=0, time_ns=0., target_ns=config['duration_ns'])
    state['pid']=os.getpid()
    atomic_json(output/'run.json', state)
    try:
        frozen = restart/'source' if resume else output/'source'
        frozen.mkdir()
        root = Path(__file__).resolve().parents[1]
        hashes = {}
        for relative in ('examples/iec_3d.py', 'plasma_bench/iec_mesh.py', 'plasma_bench/iec_field.py',
                         'plasma_bench/iec_particles.py', 'plasma_bench/iec_vacuum.py',
                         'plasma_bench/hops_drive.py', 'plasma_bench/iec_diagnostics.py',
                         'plasma_bench/iec_checkpoint.py'):
            source = root/relative
            for path in (source, source.with_suffix('.md')):
                destination = frozen/path.relative_to(root)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)
                hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
        provenance=restart/'provenance.json' if resume else output/'provenance.json'
        atomic_json(provenance, dict(hashes=hashes, numpy=np.__version__, pid=os.getpid(),
                    equations='3D electrostatic Vlasov-Poisson; linear tetra FEM; kick-drift-kick',
                    scope='collisionless source-driven development run; no neutral chemistry or fusion'))
        if resume:
            mesh=dict(np.load(output/'mesh.npz'))
        elif mesh_path:
            shutil.copy2(mesh_path, output/'mesh.npz')
            mesh = dict(np.load(output/'mesh.npz'))
            if not (np.isclose(float(mesh['size']), config['mesh_mm']*1e-3) and
                    np.isclose(float(mesh['wire_size']), config['wire_mesh_mm']*1e-3)):
                raise ValueError('Cached mesh does not match declared resolution')
        else:
            mesh = build_mesh(output/'mesh.npz', size=config['mesh_mm']*1e-3,
                              wire_size=config['wire_mesh_mm']*1e-3)
        state['status'] = 'initializing'
        atomic_json(output/'run.json', state)
        field = Field(mesh)
        rng = np.random.default_rng(config['seed'])
        p = load_plasma(field, rng, config['particles'], config['seed_radius_m'], config['density_m3'],
                        config['electron_ev'], config['ion_ev'], config['ion_amu'])
        emitter = WallEmitter(field)
        drive = CageDrive(config)
        cage_voltage = drive.initial()
        steps = int(np.ceil(config['duration_ns']*1000/config['dt_ps']))
        dt = config['duration_ns']*1e-9/steps
        source_stop_step = steps
        if config.get('source_stop_ns') is not None:
            stop_ns = config['source_stop_ns']
            if not np.isfinite(stop_ns) or not 0 <= stop_ns <= config['duration_ns']:
                raise ValueError('Source stop must be inside the physical run')
            source_stop_step = int(round(stop_ns/config['duration_ns']*steps))
            if not np.isclose(source_stop_step*config['duration_ns']/steps, stop_ns, rtol=0, atol=1e-9):
                raise ValueError('Source stop must align with an integration boundary')
        interval = max(1, int(np.ceil(steps/config['saved_intervals'])))
        charge = field.deposit(p.x, p.cells, p.charge)
        phi, electric, electrodes, energy, space = field.solve(charge, cage_voltage)
        initial = p.kinetic().sum()+energy
        loss_energy, injected_energy, source_work = 0., 0., 0.
        injected_charge = 0.
        lost_charge = np.zeros((2, 2))
        lost_count = np.zeros((2, 2), int)
        time_s = 0.
        history = []
        start_step, elapsed_offset = 0, 0.
        mesh_sha256, physics_hashes = digest(output/'mesh.npz'),source_hashes()
        if checkpoint is not None:
            meta, arrays, p = load_checkpoint(checkpoint,config,output/'mesh.npz')
            if meta['dt'] != dt or meta['scalars']['time_s'] != meta['step']*dt:
                raise ValueError('Checkpoint timestep/time mismatch')
            start_step, elapsed_offset = meta['step'],meta['elapsed_s']
            if not 0 <= start_step <= steps:
                raise ValueError('Checkpoint outside requested time interval')
            rng.bit_generator.state=meta['rng']
            vars(drive).update(meta['drive'])
            phi,charge,electric,electrodes,lost_charge,lost_count=(arrays[k] for k in
                ('phi','charge','electric','electrodes','lost_charge','lost_count'))
            s=meta['scalars']
            initial,loss_energy,injected_energy,source_work,injected_charge,time_s,energy,space,cage_voltage=(s[k] for k in
                ('initial','loss_energy','injected_energy','source_work','injected_charge','time_s','energy','space','cage_voltage'))
            history=meta['history']
            if resume and json.loads((output/'history.json').read_text()) != history:
                raise ValueError('Checkpoint must match the complete saved history; no truncation allowed')
            state.update(time_ns=time_s*1e9,step=start_step,steps=steps,frames=len(history),
                         resumed_from_us=time_s*1e6)
        # A physical midplane cut through the same tetrahedral solution.
        axis = np.linspace(-.19, .19, 45)
        gx, gz = np.meshgrid(axis, axis)
        slice_x = np.column_stack((gx.ravel(), np.zeros(gx.size), gz.ravel()))
        slice_cells = field.locate(slice_x)
        valid = slice_cells >= 0
        slice_weights = field.weights(slice_x[valid], slice_cells[valid])
        geometry = dict(vertices=mesh['vertices'].tolist(), edges=mesh['edges'].tolist(),
                        wall_radius=float(mesh['wall_radius']), wire_radius=float(mesh['wire_radius']),
                        nodes=len(field.nodes), tetrahedra=len(field.tets),
                        capacitance_pf=field.capacitance*1e12, slice_axis_m=axis.tolist())
        atomic_json(output/'geometry.json', geometry)
        state['status'] = 'running'
        def save(step):
            kinetic = p.kinetic()
            ke = float(kinetic[p.species == 0].sum())
            ki = float(kinetic[p.species == 1].sum())
            residual = ke+ki+energy+loss_energy-injected_energy-source_work-initial
            scale = max(abs(initial), abs(source_work), injected_energy, 1e-30)
            core = np.linalg.norm(p.x, axis=1) < .02
            electron_core = int(np.sum(core & (p.species == 0)))
            ion_core = int(np.sum(core & (p.species == 1)))
            centre_cell = field.locate(np.zeros((1, 3)))
            centre = float(np.dot(field.weights(np.zeros((1, 3)), centre_cell)[0], phi[field.tets[centre_cell[0]]]))
            motion = ion_motion(p)
            row = dict(time_ns=time_s*1e9, time_us=time_s*1e6, voltage_v=cage_voltage, electrons=int(np.sum(p.species == 0)),
                       ions=int(np.sum(p.species == 1)), electron_kinetic_j=ke, ion_kinetic_j=ki,
                       field_j=energy, space_charge_j=space, source_work_j=source_work,
                       injected_j=injected_energy, absorbed_j=loss_energy, energy_residual_j=float(residual),
                       energy_residual_fraction=float(residual/scale), centre_potential_v=centre,
                       minimum_potential_v=float(phi.min()), maximum_potential_v=float(phi.max()),
                       core_electrons=electron_core, core_ions=ion_core,
                       collected_charge_c=lost_charge.tolist(), collected_macroparticles=lost_count.tolist(),
                       emitted_charge_c=injected_charge, electrode_charge_c=electrodes.tolist(),
                       total_charge_residual_c=float(p.charge.sum()+lost_charge.sum()+injected_charge),
                       macro_weight=p.weight, ion_motion=motion, drive=drive.diagnostics(),
                       elapsed_s=elapsed_offset+time.monotonic()-started)
            history.append(row)
            # Display sampling uses an independent RNG and never changes the dynamics.
            display_rng = np.random.default_rng(config['seed']+step+1000000)
            selected = np.sort(display_rng.choice(len(p.x), min(len(p.x), 12000), replace=False))
            values = np.full(len(slice_x), np.nan)
            values[valid] = np.einsum('pi,pi->p', slice_weights, phi[field.tets[slice_cells[valid]]])
            atomic_json(output/f'frame_{len(history)-1:04d}.json', dict(diagnostics=row,
                positions_m=p.x[selected].tolist(), species=p.species[selected].tolist(),
                energy_ev=(kinetic[selected]/(p.weight*qe)).tolist(),
                potential_slice_v=[None if not np.isfinite(v) else float(v) for v in values]))
            snapshot=output/f'state_{len(history)-1:04d}.npz'
            if snapshot.exists():
                raise RuntimeError('Refusing to overwrite a saved scientific state')
            temporary=snapshot.with_suffix('.tmp.npz')
            np.savez_compressed(temporary, time_s=time_s,
                phi=phi, charge=charge, x=p.x, v=p.v, species=p.species,
                cells=p.cells, macro_weight=p.weight, ion_mass=p.ion_mass)
            temporary.replace(snapshot)
            state.update(time_ns=row['time_ns'], frames=len(history), step=step, steps=steps,
                         elapsed_s=row['elapsed_s'], diagnostics=row)
            atomic_json(output/'history.json', history)
            atomic_json(output/'drive_switches.json', drive.switches)
            runtime=dict(p=p,drive=drive,rng=rng,history=history,mesh_sha256=mesh_sha256,
                physics_hashes=physics_hashes,phi=phi,charge=charge,electric=electric,
                electrodes=electrodes,lost_charge=lost_charge,lost_count=lost_count,
                initial=initial,loss_energy=loss_energy,injected_energy=injected_energy,
                source_work=source_work,injected_charge=injected_charge,time_s=time_s,
                energy=energy,space=space,cage_voltage=cage_voltage)
            write_checkpoint(output/'checkpoint.npz',config,step,dt,runtime)
            atomic_json(output/'run.json', state)
            print(json.dumps({k:row[k] for k in ('time_ns', 'electrons', 'ions', 'energy_residual_fraction', 'elapsed_s')}), flush=True)
        if checkpoint is None:
            save(0)
        finish=steps if stop_step is None else min(steps,stop_step)
        if finish < start_step:
            raise ValueError('Requested stop precedes checkpoint')
        for step in range(start_step+1, finish+1):
            if (output/'STOP').exists():
                state['status'] = 'cancelled'
                if abs(history[-1]['time_ns']-time_s*1e9) > 1e-12:
                    save(step-1)
                break
            old_voltage = cage_voltage
            old_electrodes = electrodes.copy()
            # The command observes only the old state, never this step's outcome.
            if drive.mode == 'legacy':
                cage_voltage = voltage(config, step*dt)
                drive.value = cage_voltage
            else:
                radius, radial_speed, ion_count = ion_sensor(p)
                cage_voltage = drive.advance(time_s, dt, radius, radial_speed, ion_count)
            if step <= source_stop_step:
                injected, emitted = emitter.inject(p, rng, config['source_ua']*1e-6, dt, config['source_ev'])
            else:
                injected, emitted = 0., 0.
            injected_energy += injected
            injected_charge += emitted
            collection = np.array([0., emitted])  # Removing negative charge from the grounded wall.
            p.v += .5*dt*(p.charge/p.mass)[:, None]*electric[p.cells]
            target = p.x+dt*p.v
            p.x, p.cells, hit = field.drift(p.x, target, p.cells)
            if np.any(hit == 3):
                raise RuntimeError('Unresolved tetrahedral traversal; physical run stopped')
            kinetic = p.kinetic()
            for boundary in (1, 2):
                loss = hit == boundary
                collection[boundary-1] += p.charge[loss].sum()
                loss_energy += float(kinetic[loss].sum())
                for species in (0, 1):
                    selected = loss & (p.species == species)
                    lost_count[boundary-1, species] += int(selected.sum())
                    lost_charge[boundary-1, species] += p.charge[selected].sum()
            p.retain(hit == 0)
            time_s = step*dt
            new_voltage = cage_voltage
            charge = field.deposit(p.x, p.cells, p.charge)
            phi, electric, electrodes, energy, space = field.solve(charge, new_voltage)
            source_work += .5*(old_voltage+new_voltage)*(electrodes[0]-old_electrodes[0]-collection[0])
            p.v += .5*dt*(p.charge/p.mass)[:, None]*electric[p.cells]
            if not np.all(np.isfinite(p.x)) or not np.all(np.isfinite(p.v)):
                raise RuntimeError('Nonfinite physical state')
            if len(p.x) > 200000:
                raise RuntimeError('Declared 200,000 macroparticle resource bound exceeded')
            if step % interval == 0 or step == finish:
                save(step)
        else:
            state['status'] = 'completed' if finish == steps else 'checkpointed'
        np.savez_compressed(output/'final_particles.npz', time_s=time_s, x=p.x, v=p.v,
                            species=p.species, cells=p.cells, macro_weight=p.weight, ion_mass=p.ion_mass)
        np.savez_compressed(output/'final_field.npz', time_s=time_s, phi=phi, charge=charge)
        debye = np.sqrt(epsilon_0*config['electron_ev']/(config['density_m3']*qe))
        state['resolution'] = dict(initial_debye_m=debye, initial_omega_pe_dt=float(dt*np.sqrt(config['density_m3']*qe**2/(epsilon_0*electron_mass))),
                                   dt_s=dt, initial_macroparticles_per_species=config['particles'],
                                   note='Initial Debye/time scales only; evolving local density and near-wire losses are not converged.')
        atomic_json(output/'run.json', state)
        return state
    except BaseException as error:
        state.update(status='failed', error=str(error), elapsed_s=time.monotonic()-started)
        atomic_json(output/'run.json', state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    parser.add_argument('--config', type=Path)
    parser.add_argument('--mesh', type=Path)
    parser.add_argument('--checkpoint', type=Path)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--stop-step', type=int)
    args = parser.parse_args()
    config = (json.loads((Path(args.output)/'run.json').read_text())['config'] if args.resume
              else DEFAULTS | (json.loads(args.config.read_text()) if args.config else {}))
    run(args.output, config, args.mesh,checkpoint=args.checkpoint,resume=args.resume,stop_step=args.stop_step)


if __name__ == '__main__':
    main()
