"""Start the 3D hybrid model (plasma_bench/hybrid3d_device.py) from a formed FRC of the 2D formation solver, and run it.

`python -m examples.hybrid3d_frc OUT RUN t_us=10 dx=0.01 dz=0.02 ppc=8 t_end_us=10 floor=0.05 [seed=0]`

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
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE
from scipy.interpolate import RegularGridInterpolator
from plasma_bench.hybrid3d import HybridParams, curl_edges
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, DeviceHybrid, GAS
from examples.nucte_compare import load, summarize, closed_region
from examples.nucte_formation import ION_MASS
from examples.formation_display import case_from_setup, geometry
from examples.nucte_formation import coil_sections


def full_flux(run, frame_time_s):
    """Total psi on the 2D run's field grid at the saved state nearest frame_time_s; also the fluid state."""
    setup, rf, zf, frames, history = load(run)
    c = case_from_setup(setup)
    i0, r_in, fluid, grid, mapping = geometry(c)
    times = [history[k]["time_s"] for k, *_ in frames]
    pick = int(np.argmin(np.abs(np.array(times)-frame_time_s)))
    k, psi_gas, rho, p = frames[pick]
    coil = coil_sections(c)
    rr, zz = np.meshgrid(grid.r, grid.z, indexing="ij")
    unit = np.array([f[0] for f in coil.field(rr, zz)])
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
                time_s=history[k]["time_s"], mass=ION_MASS[c.gas], rows=summarize(run)[3])


def build(state, dx, dz, ppc=8, floor=.05, seed=0, dt_factor=.1):
    c, grid, fluid = state["case"], state["grid"], state["fluid"]
    dev = Device(wall_radius=c.wall_radius, tube_half_length=c.half_length,
                 coil_radius=c.wall_radius+c.exterior_cells*fluid.dr, end_gap=c.end_cells*fluid.dz)
    box = device_box(dev, dx, dz)
    types = cell_types(box, dev)
    interp = RegularGridInterpolator((grid.r, grid.z), state["psi"], method="cubic", bounds_error=False, fill_value=None)
    axis_half_bz = (16*state["psi"][1]-state["psi"][2])/(12*grid.dr**2)   # psi/r^2 on the axis (psi ~ a r^2 + b r^4)
    axis = RegularGridInterpolator((grid.z,), axis_half_bz, bounds_error=False, fill_value=None)
    lo, d = np.asarray(box.lo), box.d
    xs = [lo[a]+d[a]*np.arange(box.shape[a]) for a in range(3)]
    xh = [lo[a]+d[a]*(np.arange(box.shape[a])+.5) for a in range(3)]

    def psi_over_r2(X, Y, Z):
        # beyond the 2D grid (inside the conductor) psi is held at its edge value: no flux, B = 0 there
        rt = np.hypot(X, Y)
        r = np.minimum(rt, grid.r[-1])
        zc = np.clip(Z, grid.z[0], grid.z[-1])
        val = interp(np.column_stack((r.ravel(), zc.ravel()))).reshape(r.shape)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = val/np.maximum(rt, 1e-30)**2
        small = r < .5*grid.dr
        out[small] = axis(zc[small].reshape(-1, 1)).ravel() if small.any() else out[small]
        return out
    X, Y, Z = np.meshgrid(xh[0], xs[1], xs[2], indexing="ij")
    ax = -Y*psi_over_r2(X, Y, Z)
    X, Y, Z = np.meshgrid(xs[0], xh[1], xs[2], indexing="ij")
    ay = X*psi_over_r2(X, Y, Z)
    bf = list(curl_edges(ax, ay, np.zeros(box.shape), d))
    # ions
    mass = state["mass"]
    n_cell = state["rho"]/mass
    t_ev = np.where(n_cell > 0, state["p"]/(2*np.maximum(n_cell, 1e-30)*QE), 0.)
    fr = (fluid.rf[:-1]+fluid.rf[1:])/2
    sample = lambda arr: RegularGridInterpolator((fr, fluid.z), arr, method="linear", bounds_error=False, fill_value=None)
    f_n, f_t, f_vr, f_vz = sample(n_cell), sample(t_ev), sample(state["vr"]), sample(state["vz"])
    rng = np.random.default_rng(seed)
    gi = np.argwhere(types == GAS)
    side = int(round(ppc**(1/3)))
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
        if inside.any():
            te = float(np.sum((n_cell*t_ev)[inside])/np.sum(n_cell[inside]))
    n_peak = float(n_cell.max())
    gas_face = (types == GAS) | np.roll(types == GAS, 1, axis=2)
    b_peak = float(np.max(np.abs(bf[2][gas_face])))
    dt = dt_factor/(QE*b_peak/mass)
    params = HybridParams(ion_mass=mass, te_ev=te, n_floor=floor*n_peak)
    model = DeviceHybrid(box, dev, params, bf, pos, vel, w, dt, n_ref=n_peak)
    return model, dict(te_ev=te, n_peak=n_peak, b_peak=b_peak, dt_s=dt, particles=int(len(pos)), box=list(box.shape),
                       floor_m3=floor*n_peak, device=dict(wall=dev.wall_radius, coil=dev.coil_radius, half=dev.tube_half_length,
                                                          end_gap=dev.end_gap))


def diagnostics(model):
    """Azimuthal averages and midplane azimuthal modes."""
    box = model.box
    lo, d = np.asarray(box.lo), box.d
    xh = [lo[a]+d[a]*(np.arange(box.shape[a])+.5) for a in range(3)]
    Xf, Yf = np.meshgrid(xh[0], xh[1], indexing="ij")
    rf = np.hypot(Xf, Yf)
    dr = d[0]
    edges = np.arange(0, model.device.coil_radius+dr, dr)
    idx = np.digitize(rf, edges)-1
    nb = len(edges)-1
    # azimuthally averaged Bz(r, z) on face columns, then psi(r, z) = int Bz r dr
    bz = model.b[2]
    prof = np.zeros((nb, box.shape[2]))
    for q in range(nb):
        sel = idx == q
        if sel.any():
            prof[q] = bz[sel].mean(axis=0)
    rc = .5*(edges[:-1]+edges[1:])
    psi = np.concatenate((np.zeros((1, box.shape[2])), np.cumsum(prof*rc[:, None]*dr, axis=0)))
    rnodes = edges
    zf = lo[2]+d[2]*np.arange(box.shape[2])                  # Bz faces sit at z nodes
    gas_z = np.abs(zf) < model.device.tube_half_length
    wall_row = int(np.searchsorted(rnodes, model.device.wall_radius))
    region = closed_region(rnodes[:wall_row+1], zf[gas_z], psi[:wall_row+1, gas_z])
    # midplane density modes
    xs = [lo[a]+d[a]*np.arange(box.shape[a]) for a in range(3)]
    Xn, Yn = np.meshgrid(xs[0], xs[1], indexing="ij")
    phi = np.arctan2(Yn, Xn); rn = np.hypot(Xn, Yn)
    kz = box.shape[2]//2
    nmid = model.n[:, :, kz]
    ring = (rn > .3*model.device.wall_radius) & (rn < .9*model.device.wall_radius)
    m0 = float(np.abs(np.sum(nmid[ring])))
    modes = {f"n{m}": float(np.abs(np.sum(nmid[ring]*np.exp(-1j*m*phi[ring])))/max(m0, 1e-300)) for m in (1, 2, 3)}
    return dict(time_s=model.time, steps=model.steps, substeps=model.substeps, particles=model.particles(),
                lost=dict(model.lost), kinetic_j=model.kinetic_energy(), magnetic_j=model.magnetic_energy(),
                closed=None if region is None else bool(region["closed"]),
                trapped_flux_wb=None if region is None else region["trapped_flux_wb"],
                separatrix_mid_m=None if region is None else region["separatrix_mid_m"],
                length_m=None if region is None else region["length_m"], **modes)


def main(out, run, t_us=10., dx=.01, dz=.02, ppc=8, t_end_us=10., floor=.05, seed=0, every_us=.25):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    state = full_flux(run, t_us*1e-6)
    model, info = build(state, dx, dz, ppc, floor, seed)
    info.update(run=str(run), start_time_s=state["time_s"], dx=dx, dz=dz, ppc=ppc, floor_fraction=floor)
    (out/"setup.json").write_text(json.dumps(info, indent=1))
    history = [diagnostics(model)]
    print(json.dumps(info), flush=True)
    start = time.time()
    next_mark = every_us*1e-6
    while model.time < t_end_us*1e-6:
        model.step()
        if model.time >= next_mark-1e-15:
            row = diagnostics(model)
            row["wall_s"] = time.time()-start
            t2d = state["time_s"]+model.time
            ref = min(state["rows"], key=lambda r: abs(r["time_s"]-t2d))
            row["twod"] = dict(time_s=ref["time_s"], trapped_flux_wb=ref.get("trapped_flux_wb"),
                               separatrix_mid_m=ref.get("separatrix_mid_m"), length_m=ref.get("length_m"))
            history.append(row)
            print(json.dumps({k: v for k, v in row.items() if k != "lost"}), flush=True)
            (out/"history.json").write_text(json.dumps(history, indent=1))
            next_mark += every_us*1e-6
    return history


if __name__ == "__main__":
    args = dict(a.split("=", 1) for a in sys.argv[3:])
    main(sys.argv[1], sys.argv[2], **{k: (int(v) if k in ("ppc", "seed") else float(v)) for k, v in args.items()})
