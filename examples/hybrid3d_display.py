"""3D hybrid runs (examples/hybrid3d_frc.py) in the formation viewer, beside the 2D runs.

`python -m examples.hybrid3d_display OUT RUN[::Label] ...` writes, for each 3D run, OUT/<run>/manifest.json, one JSON
and one half-float 3D current volume per saved snapshot, and rebuilds OUT/viewer.html with every run listed in
OUT/runs.json (the 2D runs of examples/formation_display.py) followed by these. The 2D overview page is left as is.

Every frame is an actual saved snapshot (snap_NNNN.npz: B on Yee faces, node density, a random ion sample):
- **Field lines** are traced in true 3D through the node-averaged B (RK2, 5 mm steps, 3 m each way) from seeds
  on six radii x eight azimuths at the midplane. They are not rotated copies (`no_copies`): the run need not be
  axisymmetric.
- **Current** J = curl B/mu0 at the nodes, in cylindrical components (Jr, Jtheta, Jz). It is rendered as a 3D volume
  (RGBA16F texture, support = nodes touching gas).
- **Ions** are the saved random sample of macro-ions (positions), coloured by speed.
- **Separatrix (gold)** is the n = 0 (azimuthally averaged) separatrix. The frame line also gives the midplane
  n = 1 and n = 2 density amplitudes and the losses from the run history.
"""
import json
import sys
from pathlib import Path
import numpy as np
import contourpy
from scipy.constants import mu_0 as MU0, e as QE
from plasma_bench.hybrid3d import faces_to_nodes, curl_faces, edges_to_nodes
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, GAS
from examples.formation_display import case_from_setup, ring, viewer_page, GROUP_COLORS
from examples.nucte_formation import coil_sections, ION_MASS
from examples.hybrid3d_geometry import poloidal_average, density_modes


def trace(bn, lo, d, seeds, step=.005, length=3., r_max=.165, z_max=1.13):
    """Field lines through node B (trilinear), RK2, both directions; stop outside r_max or |z| > z_max."""
    shape = np.array(bn.shape[1:])

    def field(p):
        s = (p-lo)/d
        i0 = np.clip(np.floor(s).astype(int), 0, shape-2)
        f = np.clip(s-i0, 0, 1)
        out = np.zeros_like(p)
        for c in range(8):
            o = np.array([(c >> a) & 1 for a in range(3)])
            w = np.prod(np.where(o, f, 1-f), axis=1)
            idx = i0+o
            out += w[:, None]*bn[:, idx[:, 0], idx[:, 1], idx[:, 2]].T
        mag = np.linalg.norm(out, axis=1)
        return out/np.maximum(mag, 1e-30)[:, None], mag
    paths = []
    for sign in (1., -1.):
        p = seeds.copy()
        alive = np.ones(len(p), bool)
        pts, mags = [p.copy()], [field(p)[1]]
        for _ in range(int(length/step)):
            u1, _ = field(p)
            u2, _ = field(p+.5*sign*step*u1)
            q = p+sign*step*u2
            r = np.hypot(q[:, 0], q[:, 1])
            alive &= (r < r_max) & (np.abs(q[:, 2]) < z_max)
            p = np.where(alive[:, None], q, p)
            pts.append(p.copy()); mags.append(field(p)[1])
            if not alive.any():
                break
        paths.append((np.stack(pts, axis=1), np.stack(mags, axis=1)))
    out = []
    (fp, fm), (bp, bm) = paths
    for k in range(len(seeds)):
        line = np.concatenate((bp[k][::-1], fp[k][1:]))
        mag = np.concatenate((bm[k][::-1], fm[k][1:]))
        keep = np.concatenate(([True], np.any(np.diff(line, axis=0) != 0, axis=1)))
        line, mag = line[keep], mag[keep]
        if len(line) > 2:
            sel = np.unique(np.r_[np.arange(0, len(line), max(1, len(line)//150)), len(line)-1])
            out.append((line[sel], mag[sel]))
    return out


def run_display_3d(run, out, label=None):
    run = Path(run)
    info = json.loads((run/"setup.json").read_text())
    history = json.loads((run/"history.json").read_text()) if (run/"history.json").exists() else []
    dev = Device(info["device"]["wall"], info["device"]["half"], info["device"]["coil"], info["device"]["end_gap"])
    box = device_box(dev, info["dx"], info["dz"])
    types = cell_types(box, dev)
    lo, d = np.asarray(box.lo), box.d
    xs = [lo[a]+d[a]*np.arange(box.shape[a]) for a in range(3)]
    X, Y = np.meshgrid(xs[0], xs[1], indexing="ij")
    R, PHI = np.hypot(X, Y), np.arctan2(Y, X)
    gas = (types == GAS).astype(float)
    frac = gas.copy()
    for a in range(3):
        frac = .5*(frac+np.roll(frac, 1, axis=a))
    support = (frac > 0).astype(np.float32)
    setup2d = json.loads((Path(info["run"])/"setup.json").read_text())
    c = case_from_setup(setup2d)
    extent = 1.05*max(c.half_length, max(abs(s.z)+s.width/2 for s in coil_sections(c).straps))
    rr = np.array([.03, .06, .09, .115, .13, .145]); aa = np.arange(8)*np.pi/4
    seeds = np.array([[r*np.cos(a), r*np.sin(a), 0.] for r in rr for a in aa])
    snaps = sorted(run.glob("snap_*.npz"))
    if not snaps:
        raise ValueError(f"{run}: no snap_NNNN.npz (run examples/hybrid3d_frc.py with snap_us > 0)")
    (Path(out)/run.name).mkdir(parents=True, exist_ok=True)
    frames, peak_b, peak_j, vmax = [], 0., np.zeros(4), 0.
    t0 = info.get("start_time_s", 0.)
    for k, path in enumerate(snaps):
        s = np.load(path)
        b = [s["bx"].astype(float), s["by"].astype(float), s["bz"].astype(float)]
        bn = faces_to_nodes(*b)
        jn = edges_to_nodes(*[q/MU0 for q in curl_faces(*b, d)])
        cphi, sphi = np.cos(PHI)[:, :, None], np.sin(PHI)[:, :, None]
        jr = jn[0]*cphi+jn[1]*sphi; jt = -jn[0]*sphi+jn[1]*cphi; jz = jn[2]
        jmag = np.sqrt(jr**2+jt**2+jz**2)*support
        jmax = [float(jmag.max()), float(np.abs(jr*support).max()), float(np.abs(jt*support).max()), float(np.abs(jz*support).max())]
        value_scale = jmax[0] or 1.
        vol = np.stack((jr*support/value_scale, jt*support/value_scale,
                        jz*support/value_scale, support), axis=-1).astype('<f2')
        if not np.isfinite(vol).all():
            raise ValueError(f"{path}: nonfinite current volume")
        name = f"{run.name}/vol_{k:04d}.f16"
        vol.tofile(Path(out)/name)
        bmag = np.linalg.norm(bn, axis=0)
        gb = float(bmag[support > 0].max())
        lines = []
        for line, mag in trace(bn, lo, d, seeds, r_max=dev.coil_radius, z_max=dev.tube_half_length+dev.end_gap):
            pts = (np.column_stack((line[:, 2], line[:, 0], line[:, 1]))/extent).round(4).tolist()
            lines.append(dict(points=pts, field=[float(f"{m:.3g}") for m in mag]))
        rn, zn, psi, region = poloidal_average(box, dev, b[2])
        modes = density_modes(box, dev, s['n'])
        sep = []
        if region is not None:
            gen = contourpy.contour_generator(x=np.broadcast_to(zn, psi.shape),
                                              y=np.broadcast_to(rn[:, None], psi.shape), z=psi)
            o = float(psi.min())
            for pth in gen.lines(-1e-6*abs(o)):
                if len(pth) > 3:
                    for a in np.arange(4)*np.pi/2:
                        sep.append((np.column_stack((pth[:, 0], pth[:, 1]*np.cos(a), pth[:, 1]*np.sin(a)))/extent).round(4).tolist())
        ions, ti = [], None
        if "ion_x" in s.files:
            ix, iv = s["ion_x"].astype(float), s["ion_v"].astype(float)
            sp = np.linalg.norm(iv, axis=1)
            vmax = max(vmax, float(np.percentile(sp, 99)) if len(sp) else 0.)
            ions = np.column_stack((ix[:, 2]/extent, ix[:, 0]/extent, ix[:, 1]/extent, sp)).round(4).tolist()
            ti = float(np.mean(np.sum((iv-iv.mean(axis=0))**2, axis=1))*info.get("ion_mass", ION_MASS[c.gas])/(3*QE)) if len(iv) else None
        t = float(s["time_s"])
        h = next((q for q in history if abs(q["time_s"]-t) < 1e-12), {})
        exact = region or {}
        summary = ("3D hybrid · n=0: " + (f"trapped flux {1e3*exact['trapped_flux_wb']:.2f} mWb, r_s {100*exact['separatrix_mid_m']:.1f} cm, length {100*exact['length_m']:.0f} cm"
                    if exact.get('closed') else "no closed flux")
                   + f" · midplane density modes n=1 {100*modes['n1']:.1f}%, n=2 {100*modes['n2']:.1f}%"
                   + (f" · loaded ions lost {100*h['lost']['count']/max(h['particles']+h['lost']['count'], 1e-30):.1f}%" if h.get("lost") else ""))
        ijk = np.argwhere(support > 0)[::max(1, int(np.count_nonzero(support))//3000)]
        cloud = [[xs[2][z]/extent, xs[0][x]/extent, xs[1][y]/extent,
                  float(jmag[x,y,z]), float(jr[x,y,z]), float(jt[x,y,z]), float(jz[x,y,z])]
                 for x,y,z in ijk]
        frame = dict(step=int(round(t/info['dt_s'])), time_us=(t0+t)*1e6, lines=lines, current=cloud, outline=[], separatrix=sep, ions=ions,
                     bmax=gb, jmax=jmax, temperature=[info.get("te_ev"), ti] if ions else [info.get("te_ev"), None], summary=summary,
                     volume3d=dict(url=name, shape=list(box.shape), value_scale=value_scale,
                                   bounds=[xs[2][0]/extent, xs[0][0]/extent, xs[1][0]/extent, xs[2][-1]/extent, xs[0][-1]/extent, xs[1][-1]/extent]))
        if exact.get("separatrix_mid_m"):
            frame.update(midplane_radius_m=exact["separatrix_mid_m"], maximum_radius_m=float(np.nanmax(np.asarray(exact['profile'], float))))
        heavy = ("lines", "current", "outline", "separatrix", "ions")
        fname = f"{run.name}/frame_{k:04d}.json"
        (Path(out)/fname).write_text(json.dumps({q: frame[q] for q in heavy}, separators=(",", ":")))
        frames.append(dict({q: v for q, v in frame.items() if q not in heavy}, url=fname))
        peak_b = max(peak_b, gb); peak_j = np.maximum(peak_j, jmax)
        print(f"{run.name}: frame {k+1}/{len(snaps)} t={t*1e6:.2f} us", flush=True)
    chamber = [ring(dev.wall_radius, z, extent) for z in (-dev.tube_half_length, 0., dev.tube_half_length)]
    chamber += [(np.array([[-dev.tube_half_length, dev.wall_radius*np.cos(a), dev.wall_radius*np.sin(a)],
                           [dev.tube_half_length, dev.wall_radius*np.cos(a), dev.wall_radius*np.sin(a)]])/extent).round(5).tolist()
                for a in np.arange(8)*np.pi/4]
    coil = coil_sections(c)
    coil_groups = [dict(label="Main coil (perfect conductor in 3D)", color=GROUP_COLORS[0],
                        paths=[ring(st.radius, st.z, extent, 61) for st in coil.straps if st.section == 0])]
    title = label or run.name
    description = (f"3D hybrid run (kinetic ions, fluid electrons with the Hall term) started at {t0*1e6:.1f} µs from the 2D run "
                   f"{Path(info['run']).name}; grid {info['box'][0]}×{info['box'][1]}×{info['box'][2]} ({1e3*info['dx']:.0f} mm × {1e3*info['dz']:.0f} mm), "
                   f"{info['particles']:,} macro-ions, isothermal electrons at {info['te_ev']:.0f} eV")
    scope = (description+f". Declared density floor in Ohm's law {info['floor_m3']:.2g} m⁻³ ({100*info['floor_fraction']:.0f}% of peak). "
             "Field lines are traced in 3D (not rotated copies); current is the full 3D J = curl B/μ0 in cylindrical components; "
             "ions are independent random samples, not persistent particle trails. Losses count the entire loaded gas, including ions outside the initial separatrix. "
             "The ion velocity-spread equivalent includes bulk-flow variation and is unweighted. "
             "The gold separatrix is the n = 0 (azimuthal average). "
             "M3 of docs/HYBRID_3D.md: the 3D model is newly verified (waves, vacuum gap, compressional eigenmode); results are preliminary.")
    data = dict(frames=frames, detailed=True, no_copies=True, bmax=peak_b, jmax=peak_j.tolist(), ion_speed_max=vmax or 1.,
                initial_frame=0, playback_fps=10, initial_yaw=-.35, initial_pitch=.28, initial_zoom=3.,
                chamber=chamber, coils=[], coil_groups=coil_groups, physical_scale_m=extent,
                status_note=f"{title} · saved 3D hybrid states", scope=scope, temperature_digits=3,
                drive_label="Coil", motion_label="Separatrix radius (n=0) — midplane", title=title,
                description=description+". Gold: n = 0 separatrix. Dots (Ions checkbox): sampled macro-ions.")
    data.update(solver_dt_s=info['dt_s'], solver_steps=frames[-1]['step'],
                model_label="3D hybrid · kinetic ions", temperature_label="Te / sampled ion velocity spread",
                model_notice="Preliminary 3D continuation of a formed FRC. Stability and confinement have not been established.",
                rendering_note="Current is sampled from the saved 3D field, without rotational copies. Colors show A/m²; opacity is illustrative, not emitted light. Ion dots are independent samples, not trajectories.")
    (Path(out)/run.name/"manifest.json").write_text(json.dumps(data, separators=(",", ":")))
    return dict(run=run.name, label=title, description=description, manifest=f"{run.name}/manifest.json", frames=len(frames),
                inner_radius_m=0., rows=[], flux_efold_s=None, fit_span_s=0., hybrid3d=True)


def main(out, runs):
    out = Path(out)
    existing = json.loads((out/"runs.json").read_text()) if (out/"runs.json").exists() else []
    replaced = {Path(spec.partition('::')[0]).name for spec in runs}
    existing = [s for s in existing if s['run'] not in replaced]
    added = []
    for spec in runs:
        run, _, label = spec.partition("::")
        added.append(run_display_3d(run, out, label or None))
    allruns = existing+added
    (out/"runs.json").write_text(json.dumps(allruns, indent=1, default=float))
    viewer_page(out, allruns)
    html = (out/"viewer.html").read_text()
    if 'id="ions"' not in html:
        html = html.replace('<label><input type="checkbox" id="separatrix" checked> Separatrix</label>',
                            '<label><input type="checkbox" id="separatrix" checked> Separatrix</label>'
                            '<label><input type="checkbox" id="ions" checked> Ions (3D runs)</label>')
    if 'id="ion-legend"' not in html:
        html = html.replace('<p id="frame-summary"></p>', '<p id="frame-summary"></p><p id="ion-legend" class="quiet"></p>')
    (out/"viewer.html").write_text(html)
    return added


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
