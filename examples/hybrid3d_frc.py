"""Start the 3D hybrid model (plasma_bench/hybrid3d_device.py) from a formed FRC of the 2D formation solver, and run it.

`python -m examples.hybrid3d_frc OUT RUN t_us=10 dx=0.01 dz=0.02 ppc=8 t_end_us=10 floor=0.05 [seed=0 snap_us=0.1]`
(snap_us: interval of OUT/snap_NNNN.npz field + ion-sample snapshots for examples/hybrid3d_display.py; 0 disables)

Loading (M3 of docs/HYBRID_3D.md):
- **Field.** The 2D run's total poloidal flux at the chosen saved state is rebuilt on its full field grid (gas plus
  vacuum), as examples/formation_display.py does: the saved gas-node flux, minus the section coil flux of the
  saved vacuum currents, extended by the run's own exterior solve (CouplingMesh.harmonic_extension), plus that
  coil flux. The 3D B is the Yee curl of A_theta = psi/r placed on the edges, so div B = 0 exactly. Near the
  axis psi/r^2 is used, which is regular.
- **Ions.** `ppc` macro-ions on a lattice in every gas cell, weighted by the local density. Their velocity is
  the 2D fluid velocity (v_r, v_z) plus a Maxwellian at T_i = p/(2 n e) (the 2D model has Te = Ti).
- **Electrons.** Isothermal at the density-weighted mean T inside the separatrix.
- **Declared density floor** in Ohm's law: `floor` times the peak density.
- **Mapping.** The 2D boundary (fixed induced flux at the coil radius and two cells beyond the ends) maps onto
  the 3D conductor at the coil and the end planes. Exact vacuum fills the gap in both.

Diagnostics every `every_us`:
- the azimuthal average of B_z, giving psi(r, z), the separatrix and the trapped flux;
- the n = 1 and n = 2 azimuthal Fourier amplitudes of density at the midplane, relative to n = 0;
- particles and losses, ion kinetic and magnetic energy.
The 2D run's own trapped flux and separatrix at the same times are listed for comparison.
"""
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE
from scipy.interpolate import RegularGridInterpolator
from plasma_bench.hybrid3d import HybridParams, curl_edges
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, DeviceHybrid, GAS
from plasma_bench.hybrid3d_energy_electrons import EnergyElectrons
from examples.nucte_compare import load, summarize
from examples.nucte_formation import ION_MASS
from examples.formation_display import case_from_setup, geometry
from examples.nucte_formation import coil_sections
from examples.hybrid3d_geometry import poloidal_average, density_modes
from examples.hybrid3d_sources import flux_edges, BankWaveform
from examples.hybrid3d_circuit import CircuitDrive


def full_flux(run, frame_time_s):
    """Total psi on the 2D run's field grid at the saved state nearest frame_time_s; also the fluid state."""
    setup, rf, zf, frames, history = load(run)
    c = case_from_setup(setup)
    if (Path(run)/'toroidal.npz').exists():
        raise ValueError('3D initialization currently maps poloidal sources only; toroidal field/rotation cannot be discarded')
    i0, r_in, fluid, grid, mapping = geometry(c)
    coil = coil_sections(c)
    rr, zz = np.meshgrid(grid.r, grid.z, indexing='ij')
    unit = np.array([f[0] for f in coil.field(rr, zz)])
    if frame_time_s == 0:
        # Old archives begin at 0.1 us. Rebuild the declared initial fill,
        # rather than relabeling that evolved state as time zero.
        if history[0]['time_s'] != 0:
            raise ValueError('Initial coil-current state at t=0 is unavailable')
        fraction = np.ones((c.nr, c.nz))
        if c.pi_half_length > 0:
            shape = .5*(1-np.tanh((abs(fluid.z)-c.pi_half_length)/c.pi_edge_width))
            residual = max(c.floor_fraction, 1e-3)
            fraction *= residual+(1-residual)*shape[None, :]
        rho = setup['density_kg_m3']*fraction
        p = setup['pressure_pa']*fraction
        return dict(case=c, setup=setup, grid=grid,
                    psi=np.einsum('s,srz->rz', history[0]['circuit']['vacuum_current_a'], unit),
                    fluid=fluid, rho=rho, p=p, vr=np.zeros_like(rho), vz=np.zeros_like(rho),
                    time_s=0., mass=ION_MASS[c.gas], rows=summarize(run)[3], unit_flux=unit,
                    drive_history=history, inner_radius=r_in)
    times = [history[k]["time_s"] for k, *_ in frames]
    pick = int(np.argmin(np.abs(np.array(times)-frame_time_s)))
    k, psi_gas, rho, p = frames[pick]
    vac = np.asarray(history[k]["circuit"]["vacuum_current_a"])
    source = np.tensordot(vac, unit, axes=1)
    induced = np.zeros(grid.shape)
    rows, cols = slice(i0, i0+c.nr+1), slice(mapping.z_offset, mapping.z_offset+c.nz+1)
    induced[rows, cols] = psi_gas-source[rows, cols]
    psi = mapping.harmonic_extension(induced)+source
    vel = np.load(Path(run)/"velocity.npz")
    key = f"t{k:04d}" if f"t{k:04d}" in vel.files else f"t{k:02d}"
    vv = vel[key]
    vr, vz = vv[:c.nr*c.nz].reshape(c.nr, c.nz), vv[c.nr*c.nz:].reshape(c.nr, c.nz)
    return dict(case=c, setup=setup, grid=grid, psi=psi, fluid=fluid, rho=rho, p=p, vr=vr, vz=vz,
                time_s=history[k]["time_s"], mass=ION_MASS[c.gas], rows=summarize(run)[3],
                unit_flux=unit, drive_history=history, inner_radius=r_in)


def build(state, dx, dz, ppc=8, floor=.05, seed=0, dt_factor=.1, drive=0, duration_s=40e-6, gamma_e=0., dt_ns=0., thermal=0, gpu=0, fused_arrays=False):
    if not np.isfinite(ppc) or ppc < 1 or ppc != int(ppc):
        raise ValueError('ppc must be a positive perfect cube: 1, 8, 27, 64, ...')
    side = int(round(ppc**(1/3)))
    if side**3 != ppc:
        raise ValueError(f'ppc={ppc} is not a perfect cube; the cubic loader would use {side**3}. '
                         'Choose an exact count such as 8, 27 or 64.')
    c, grid, fluid = state["case"], state["grid"], state["fluid"]
    dev = Device(wall_radius=c.wall_radius, tube_half_length=c.half_length,
                 coil_radius=c.wall_radius+c.exterior_cells*fluid.dr, end_gap=c.end_cells*fluid.dz,
                 inner_radius=state['inner_radius'])
    box = device_box(dev, dx, dz)
    types = cell_types(box, dev)
    lo, d = np.asarray(box.lo), box.d
    bf = list(curl_edges(*flux_edges(box, grid, state['psi']), d))
    source = (CircuitDrive(box, state, duration_s, types == GAS) if drive == 2 else
              BankWaveform(box, grid, state['unit_flux'], state['drive_history'], state['time_s'], types == GAS) if drive else None)
    if source is not None:
        plasma_psi = state['psi']-np.einsum('s,srz->rz', source.coefficients(0.), state['unit_flux'])
        bf = (list(source(0.)) if state['time_s'] == 0 else
              list(np.asarray(curl_edges(*flux_edges(box, grid, plasma_psi), d))+source(0.)))
    if c.aux_coils and not drive:
        raise ValueError('Auxiliary coils require drive=1 so vacuum projection retains their currents')
    # ions
    mass = state["mass"]
    n_cell = state["rho"]/mass
    t_ev = np.where(n_cell > 0, state["p"]/(2*np.maximum(n_cell, 1e-30)*QE), 0.)
    fr = (fluid.rf[:-1]+fluid.rf[1:])/2
    sample = lambda arr: RegularGridInterpolator((fr, fluid.z), arr, method="linear", bounds_error=False, fill_value=None)
    f_n, f_t, f_vr, f_vz = sample(n_cell), sample(t_ev), sample(state["vr"]), sample(state["vz"])
    rng = np.random.default_rng(seed)
    gi = np.argwhere(types == GAS)
    offs = np.array([((i+.5)/side, (j+.5)/side, (l+.5)/side) for i in range(side) for j in range(side) for l in range(side)])
    pos = (lo+(gi[:, None, :]+offs[None, :, :])*d).reshape(-1, 3)
    r = np.hypot(pos[:, 0], pos[:, 1])
    rz = np.column_stack((np.clip(r, fr[0], fr[-1]), np.clip(pos[:, 2], fluid.z[0], fluid.z[-1])))
    n_p = np.maximum(f_n(rz), 0.); t_p = np.maximum(f_t(rz), .5)
    vr, vz = f_vr(rz), f_vz(rz)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx, sx = np.where(r > 0, pos[:, 0]/r, 1.), np.where(r > 0, pos[:, 1]/r, 0.)
    vel = np.column_stack((vr*cx, vr*sx, vz))+rng.normal(size=pos.shape)*np.sqrt(t_p*QE/mass)[:, None]
    w = n_p*box.cell_volume/len(offs)
    keep = w > 0
    pos, vel, w = pos[keep], vel[keep], w[keep]
    # electron temperature: density-weighted mean T inside the 2D separatrix at this time
    ref = min(state["rows"], key=lambda q: abs(q["time_s"]-state["time_s"]))
    te = float(np.sum(n_cell*t_ev)/np.sum(n_cell))
    if ref.get("profile") is not None:
        prof = np.array([np.nan if v is None else v for v in ref["profile"]], float)
        inside = np.zeros_like(n_cell, bool)
        for j in range(max(ref["lo"], 0), min(ref["hi"], len(prof)-1)):
            rs = np.nanmin([prof[j], prof[j+1]])
            inside[:, j] = fr < rs
            if ref.get('inner_profile') is not None:
                inner = np.asarray(ref['inner_profile'], float)
                inside[:, j] &= fr > np.nanmax([inner[j], inner[j+1]])
        if inside.any():
            te = float(np.sum((n_cell*t_ev)[inside])/np.sum(n_cell[inside]))
    n_peak = float(n_cell.max())
    gas_face = (types == GAS) | np.roll(types == GAS, 1, axis=2)
    b_peak = float(np.max(np.abs(bf[2][gas_face])))
    dt = dt_factor/(QE*b_peak/mass)
    if dt_ns > 0:
        dt = min(dt, dt_ns*1e-9)
    gamma = gamma_e if gamma_e > 1 else 5/3 if thermal else None
    params = HybridParams(ion_mass=mass, te_ev=te, gamma_e=gamma, n_floor=floor*n_peak)
    electron_options = dict(closure='energy', gamma=gamma, resistivity='spitzer', coulomb_log=10., transport_cfl=.8) if thermal else None
    if fused_arrays:
        if not thermal:
            raise ValueError('Fused energy arrays require thermal=1')
        electron_options['fused_arrays'] = True
    if gpu:
        if not thermal:
            raise ValueError('GPU prototype requires thermal=1')
        electron_options['backend'] = 'metal'
        electron_options['gpu_moments'] = gpu >= 2
    model = DeviceHybrid(box, dev, params, bf, pos, vel, w, dt, n_ref=n_peak, source=source, thermal=electron_options)
    return model, dict(te_ev=te, n_peak=n_peak, b_peak=b_peak, dt_s=dt, particles=int(len(pos)), box=list(box.shape),
                       particle_loading=dict(positions='cell-centred tensor-product lattice', lattice_shape=[side]*3,
                                             samples_per_gas_cell=len(offs), gas_cells=len(gi),
                                             zero_weight_samples_removed=int(len(keep)-np.count_nonzero(keep))),
                       floor_m3=floor*n_peak, device=dict(wall=dev.wall_radius, coil=dev.coil_radius, half=dev.tube_half_length,
                                                          end_gap=dev.end_gap, inner=dev.inner_radius),
                       drive='coupled_3d_banks' if drive == 2 else 'prescribed_saved_bank_currents' if drive else 'fixed_conductor_flux',
                       thermal_electrons=electron_options,
                       compute_backend=('metal_transport_and_cic_cpu64_integration' if gpu >= 2 else
                                        'metal_transport_cpu64_integration' if gpu else 'cpu_float64'),
                       electron_boundary_energy=('outgoing fluid enthalpy plus energy left on emptied nodes' if thermal else None),
                       plasma_current_support='gas_edges_excluding_conductors',
                       current_node_volume='gas_fraction',
                       symmetry='full Cartesian domain; no azimuthal averaging in evolution; axisymmetric initial fluid state')


def diagnostics(model):
    """Azimuthal averages and midplane azimuthal modes."""
    *_, region = poloidal_average(model.box, model.device, model.b[2])
    modes = density_modes(model.box, model.device, model.n)
    row = dict(time_s=model.time, steps=model.steps, substeps=model.substeps, particles=model.particles(),
                lost=dict(model.lost), kinetic_j=model.kinetic_energy(), magnetic_j=model.magnetic_energy(),
                closed=None if region is None else bool(region["closed"]),
                trapped_flux_wb=None if region is None else region["trapped_flux_wb"],
                separatrix_mid_m=None if region is None else region["separatrix_mid_m"],
                length_m=None if region is None else region["length_m"], **modes)
    pressure = model._pressure(model.n)
    number = float(np.sum(model.n*model.node_gas_fraction)*model.box.cell_volume)
    row['mean_te_ev'] = float(np.sum(pressure*model.node_gas_fraction)*model.box.cell_volume/max(number*QE, 1e-300))
    row['electron_internal_j'] = (float(np.sum(pressure*model.node_gas_fraction)*model.box.cell_volume/(model.p.gamma_e-1))
                                  if model.p.gamma_e is not None else None)
    row['max_div_B_relative'] = model.max_divergence()
    bmax = np.sqrt(sum(float(np.max(np.abs(b[active])))**2 for b, active in zip(model.b, model.face_active)))
    row['omega_dt_bound'] = float(model.qm*bmax*model.dt)
    row['cell_travel_bound'] = float(np.max(np.abs(model.v)/model.d)*model.dt) if model.v.size else 0.
    if isinstance(model.source, CircuitDrive):
        row['circuit'] = model.source.snapshot()
    if model.electrons is not None:
        row.update(ohmic_heat_j=model.electrons.ohmic_heat_j, cutoff_removed_heat_j=model.electrons.cutoff_removed_heat_j,
                   electron_transport_fraction=model.electrons.max_transport_fraction,
                   absorbed_electron_energy_j=model.electrons.absorbed_energy_j)
        if isinstance(model.electrons, EnergyElectrons):
            row.update(electron_boundary_enthalpy_j=model.electrons.enthalpy_outflow_j,
                       electron_empty_node_energy_j=model.electrons.vacuum_removed_energy_j,
                       electron_grid_balance_residual_j=model.electrons.energy_balance_residual_j,
                       electron_pressure_centering_residual_j=model.electrons.pressure_centering_residual_j,
                       electron_roundoff_energy_j=model.electrons.roundoff_energy_j)
            row.update(electron_energy_substeps=model.electrons.energy_substeps,
                       electron_energy_max_subcycles=model.electrons.max_energy_subcycles,
                       electron_energy_max_loss_fraction=model.electrons.max_energy_loss_fraction)
            if hasattr(model.electrons, 'metal'):
                row.update(metal_kernel_calls=model.electrons.metal.calls,
                           metal_synchronized_wall_s=model.electrons.metal.wall_s)
        else:
            row.update(electron_remap_energy_j=model.electrons.remap_energy_j,
                       electron_heat_deposition_residual_j=model.electrons.heat_deposition_residual_j,
                       electron_entropy_clip_count=model.electrons.entropy_clip_count)
    if isinstance(model.source, CircuitDrive) and isinstance(model.electrons, EnergyElectrons):
        response = model.response_magnetic_energy()
        accounted = (response+model.kinetic_energy()+model.electrons.energy(model.n)
                     +model.lost['energy_j']+model.electrons.absorbed_energy_j+model.electrons.cutoff_removed_heat_j)
        if not hasattr(model, '_initial_load_energy'):
            model._initial_load_energy = accounted
        row.update(response_magnetic_j=response,
                   load_energy_work_j=row['circuit'].get('source_work_j', row['circuit']['image_work_j']),
                   load_energy_residual_j=accounted-model._initial_load_energy-row['circuit'].get('source_work_j', row['circuit']['image_work_j']))
    return row


def save_snapshot(model, out, index, rng, ions=4000):
    """Fields and an ion sample for the 3D viewer (examples/hybrid3d_display.py): B on faces and node density as
    float32, and `ions` randomly chosen macro-ions (position, velocity), a point sample with no identity."""
    pick = rng.choice(model.x.shape[0], size=min(ions, model.x.shape[0]), replace=False)
    path = Path(out)/f"snap_{index:04d}.npz"
    temporary = path.with_suffix('.npz.tmp')
    source_data = {}
    if model.source is not None:
        source_data = dict(b_response=(np.asarray(model.b)-model.source(model.time)).astype(np.float32),
                           source_current_a=model.source.coefficients(model.time))
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, time_s=model.time,
                            bx=model.b[0].astype(np.float32), by=model.b[1].astype(np.float32), bz=model.b[2].astype(np.float32),
                            n=model.n.astype(np.float32), pe=model._pressure(model.n).astype(np.float32),
                            ion_x=model.x[pick].astype(np.float32), ion_v=model.v[pick].astype(np.float32), **source_data)
    temporary.replace(path)


def save_history(history, out):
    temporary = Path(out)/'history.json.tmp'
    temporary.write_text(json.dumps(history, indent=1))
    temporary.replace(Path(out)/'history.json')


def main(out, run, t_us=10., dx=.01, dz=.02, ppc=8, t_end_us=10., floor=.05, seed=0, every_us=.25, snap_us=.1, drive=0, gamma_e=0., dt_ns=0., thermal=0, gpu=0, arrays=0):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    state = full_flux(run, t_us*1e-6)
    model, info = build(state, dx, dz, ppc, floor, seed, drive=int(drive), duration_s=t_end_us*1e-6,
                        gamma_e=gamma_e, dt_ns=dt_ns, thermal=bool(thermal), gpu=int(gpu), fused_arrays=bool(arrays))
    if int(drive) == 1:
        model.source.coefficients(t_end_us*1e-6+model.dt)
        (out/'drive.json').write_text(json.dumps(dict(time_s=model.source.times.tolist(),
                                                     vacuum_current_a=model.source.currents.tolist()), indent=1))
    info.update(run=str(run), start_time_s=state["time_s"], dx=dx, dz=dz, ppc=ppc, floor_fraction=floor,
                ion_mass=state['mass'], seed=seed, requested_duration_s=t_end_us*1e-6,
                snapshot_interval_s=snap_us*1e-6, diagnostic_interval_s=every_us*1e-6,
                hybrid_parameters=asdict(model.p), cell_spacing_m=model.box.d.tolist())
    if int(drive) == 2:
        info['bank_calibration'] = model.source.calibration
    if model.source is not None:
        np.savez_compressed(out/'coil_basis.npz', b_unit=model.source.faces)
        info['current_definition'] = 'plasma current: curl(B_total-B_known_coils)/mu0'
    (out/"setup.json").write_text(json.dumps(info, indent=1))
    history = [diagnostics(model)]
    save_history(history, out)
    print(json.dumps(info), flush=True)
    start = time.time()
    next_mark = every_us*1e-6
    snap_rng = np.random.default_rng(seed+1)
    snaps, next_snap = 0, 0.
    last_snap_time = None
    if snap_us > 0:
        save_snapshot(model, out, 0, snap_rng); snaps, next_snap = 1, snap_us*1e-6
        last_snap_time = model.time
    while model.time < t_end_us*1e-6:
        model.step()
        if isinstance(model.source, CircuitDrive):
            model.source.advance(model)
        if snap_us > 0 and model.time >= next_snap-1e-15:
            save_snapshot(model, out, snaps, snap_rng); snaps += 1; next_snap += snap_us*1e-6
            last_snap_time = model.time
        if model.time >= next_mark-1e-15:
            row = diagnostics(model)
            row["wall_s"] = time.time()-start
            t2d = state["time_s"]+model.time
            ref = min(state["rows"], key=lambda r: abs(r["time_s"]-t2d))
            row["twod"] = dict(time_s=ref["time_s"], trapped_flux_wb=ref.get("trapped_flux_wb"),
                               separatrix_mid_m=ref.get("separatrix_mid_m"), length_m=ref.get("length_m"))
            history.append(row)
            print(json.dumps({k: v for k, v in row.items() if k != "lost"}), flush=True)
            save_history(history, out)
            next_mark += every_us*1e-6
    if snap_us > 0 and last_snap_time != model.time:
        save_snapshot(model, out, snaps, snap_rng)
    if history[-1]['time_s'] != model.time:
        row = diagnostics(model)
        ref = min(state['rows'], key=lambda r: abs(r['time_s']-state['time_s']-model.time))
        row['twod'] = {key: ref.get(key) for key in ('time_s', 'trapped_flux_wb', 'separatrix_mid_m', 'length_m')}
        row['wall_s'] = time.time()-start
        history.append(row)
    save_history(history, out)
    return history


if __name__ == "__main__":
    args = dict(a.split("=", 1) for a in sys.argv[3:])
    main(sys.argv[1], sys.argv[2], **{k: (int(v) if k in ("ppc", "seed") else float(v)) for k, v in args.items()})
