"""Reproducible actual ideal electron/ion flow measurements."""
from pathlib import Path
from time import perf_counter
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd import MHDParameters
from plasma_bench.toroidal_state import ToroidalState
from plasma_bench.warm_mhd import AdiabaticSpeciesMHD
from .warm_reference import WarmReference
from .warm_cases import protocol
from .warm_records import provenance, save
from .verification_records import write_record


def history_columns(rows):
    """Preserve every accepted value without repeating field labels per step."""
    return {key: history_columns([row[key] for row in rows]) if isinstance(value, dict) else [row[key] for row in rows]
            for key, value in rows[0].items()}


def initialize(case):
    ref = WarmReference(**case['reference'])
    mesh = MHDMesh(case['radius'], ref.length/2, case['nr'], case['nz'])
    means = ref.cell_means(mesh)
    values = np.zeros((mesh.nr, mesh.nz, 6))
    values[..., 0], values[..., 2], values[..., 3] = means['density'], means['momentum'], means['energy']
    r, _ = mesh.corners
    state = ToroidalState(values, ref.bias*r*r/2, 'compatible')
    pressure = np.broadcast_to(ref.fields(mesh.zf[:-1], 0)['electron_pressure'], (mesh.nr+1, mesh.nz))
    model = AdiabaticSpeciesMHD(mesh, state, pressure, ref.ion_mass,
                               MHDParameters(gamma=ref.gamma, cfl=case['cfl'], spatial_order=case['order']),
                               density_scale=case['density_scale'])
    return ref, mesh, model


def measure(case):
    ref, mesh, model = initialize(case)
    initial_state, initial_u = model.state.copy(), model.electron_entropy.copy()
    initial = model.diagnostics()
    history = [initial]
    end = case['periods']*ref.period
    steps = case['fixed_steps']
    if steps is not None:
        for _ in range(steps):
            model.step(end/steps)
            history.append(model.diagnostics())
    else:
        while model.time < end:
            remaining = end-model.time
            if remaining <= 4*np.finfo(float).eps*end:
                break
            model.step(min(remaining, model.timestep_limit()))
            history.append(model.diagnostics())
    final = model.diagnostics()
    pe, pi = model.thermal.validate(model.state, model.electron_entropy)
    true_node = ref.fields(mesh.zf[:-1], model.time)['electron_pressure']
    true_cell = ref.cell_means(mesh, model.time)
    pe_scale = np.sqrt(np.mean((true_node-ref.pe0)**2))
    error = np.sqrt(np.sum(model.thermal.mapping.volumes*(pe-true_node)**2)/np.sum(mesh.volumes))
    metrics = dict(electron_node_error_pa=float(error), electron_node_perturbation_relative_L2=float(error/pe_scale),
                   ion_cell_error_pa=float(np.sqrt(np.sum(mesh.volumes*(pi-true_cell['ion_pressure'])**2)/np.sum(mesh.volumes))),
                   density_cell_relative_L2=float(np.sqrt(np.sum(mesh.volumes*(model.state.conserved[..., 0]-true_cell['density'])**2)/np.sum(mesh.volumes))/ref.rho0),
                   mass_relative_change=(final['mass_kg']-initial['mass_kg'])/initial['mass_kg'],
                   total_energy_relative_change=(final['total_energy_j']-initial['total_energy_j'])/initial['total_energy_j'],
                   electron_log_entropy_change_per_particle=(final['electron_log_entropy_particles']-initial['electron_log_entropy_particles'])/(initial['mass_kg']/ref.ion_mass))
    arrays = dict(initial_conserved=initial_state.conserved, final_conserved=model.state.conserved,
                  initial_psi=initial_state.psi, final_psi=model.state.psi,
                  initial_entropy=initial_u, final_entropy=model.electron_entropy,
                  final_electron_pressure=pe, final_ion_pressure=pi,
                  primal_volumes=mesh.volumes, dual_volumes=model.thermal.mapping.volumes,
                  radial_faces=mesh.rf, axial_faces=mesh.zf)
    return dict(actual_time_s=model.time, end_time_s=end, accepted_steps=model.steps, initial=initial,
                final=final, history_columns=history_columns(history), metrics=metrics), arrays


def run(directory='outputs/warm_adiabatic_verification'):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source, records = provenance(), []
    for case in protocol():
        start = perf_counter()
        record, arrays = measure(case)
        records.append(save(directory, case, source, record | dict(elapsed_seconds=perf_counter()-start), arrays))
        print(case['name'], record['metrics'], flush=True)
    if provenance() != source:
        raise RuntimeError('Sources changed during warm protocol')
    write_record(directory/'measurements.json', dict(completed=True, provenance=source, cases=records))


if __name__ == '__main__':
    run()
