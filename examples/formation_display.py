"""Interactive 3D display of formation runs (examples/nucte_formation.py) in the saved-field viewer.

`python -m examples.formation_display OUT RUN[::Label] ...` writes OUT/viewer.html, one 3D viewer for all
runs with a run selector (visualizations/material_fields volume.html + viewer.js + volume.js), the per-run
manifests and saved-state frames under OUT/<run>/, and OUT/index.html: trapped flux, separatrix length and
radius against time for all runs (at most six), a table and links into the viewer. Serve OUT over HTTP.

Read-only: every frame is a saved 1 µs snapshot; nothing is evolved or interpolated in time.
- **Full poloidal flux.** The snapshot keeps the total flux on the gas nodes only. The induced part,
  total − Σ_s I_vac,s ψ_s, is extended into the vacuum (outside the wall, inside an inner tube, beyond
  the ends) by the run's own exterior solve, K u = 0 (CouplingMesh.harmonic_extension). The coil flux is
  then added back with the same saved section currents, so the gas-node flux is reproduced exactly.
- **Field lines** are ψ contours (equal flux between open lines; closed surfaces inside the separatrix)
  in four meridional planes, coloured by |B| from the same flux.
- **Current** is J_θ = −Δ*ψ/(μ0 r) from second differences of that flux on the gas nodes.
  J_r = J_z = 0: the model has no toroidal field.
- **The separatrix** (gold) and FRC quantities come from examples.nucte_compare (closed_region, or
  annular_closed_region for an inner tube).
- **Particles** (optional): RUN/particles.json from examples/particle_tracking.py, actual integrated ion
  orbits in the same saved fields, drawn as trails in the field canvas behind a "Particles" checkbox.
"""
import base64
import json
import sys
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
import contourpy
from scipy.constants import mu_0 as MU0
from scipy.interpolate import RegularGridInterpolator
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.geometry import Grid
from plasma_bench.coupling_mesh import CouplingMesh
from examples.nucte_formation import Nucte, AuxCoil, Bank, coil_sections, inner_rows, _from_dict
from examples.nucte_compare import load, summarize, decay
from examples.design_scan_summary import formed

ASSETS = Path(__file__).resolve().parents[1]/"visualizations/material_fields"
GROUP_COLORS = ("#db9160", "#96b5ef", "#c6a5d8", "#8fd19e", "#e8c46b")   # main coil first, then each AuxCoil
AZIMUTHS = np.arange(4)*np.pi/2


def case_from_setup(setup):
    case = dict(setup["case"])
    for key, kind in (("aux_coils", AuxCoil), ("extra_main_banks", Bank)):
        case[key] = tuple(_from_dict(kind, d) for d in case.get(key, ()))
    known = set(Nucte.__dataclass_fields__)
    return replace(Nucte(), **{k: v for k, v in case.items() if k in known})


def geometry(c):
    i0, r_in = inner_rows(c)
    fluid = MHDMesh(c.wall_radius, c.half_length, c.nr, c.nz, r_in)
    grid = Grid(c.wall_radius+c.exterior_cells*fluid.dr, c.half_length+c.end_cells*fluid.dz,
                i0+c.nr+c.exterior_cells+1, c.nz+2*c.end_cells+1)
    return i0, r_in, fluid, grid, CouplingMesh(fluid, grid, boundary="fixed")


def ring(radius, z, extent, n=97):
    a = np.linspace(0, 2*np.pi, n)
    return (np.column_stack((np.full(n, z), radius*np.cos(a), radius*np.sin(a)))/extent).round(5).tolist()


def rotate(path_zr, extent, azimuths=AZIMUTHS):
    z, r = path_zr[:, 0], path_zr[:, 1]
    return [(np.column_stack((z, r*np.cos(a), r*np.sin(a)))/extent).round(4).tolist() for a in azimuths]


def decimate(path, limit=60):
    if len(path) <= limit:
        return path
    keep = np.unique(np.r_[np.arange(0, len(path), int(np.ceil(len(path)/limit))), len(path)-1])
    return path[keep]


def particle_trails(run, out, extent):
    """Copy RUN/particles.json (examples/particle_tracking.py: actual integrated orbits, metres) to
    OUT/<run>/particles.json in viewer coordinates [z, x, y]/extent; return its relative url or None."""
    source = Path(run)/"particles.json"
    if not source.exists():
        return None
    data = json.loads(source.read_text())
    for trail_set in data["sets"]:
        for trail in trail_set["trails"]:
            p = np.asarray(trail["points"], float)
            trail["points"] = (np.column_stack((p[:, 2], p[:, 0], p[:, 1]))/extent).round(4).tolist()
    name = f"{Path(run).name}/particles.json"
    (Path(out)/name).write_text(json.dumps(data, separators=(",", ":")))
    return name


def fields(grid, psi):
    """Bz = (1/r) dpsi/dr, Br = -(1/r) dpsi/dz on the full grid; regular axis Bz = 2a for psi = a r^2 + b r^4."""
    r = grid.r[:, None]
    dpsi_dr = np.gradient(psi, grid.dr, axis=0, edge_order=2)
    dpsi_dz = np.gradient(psi, grid.dz, axis=1, edge_order=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        bz, br = dpsi_dr/r, -dpsi_dz/r
    bz[0] = (16*psi[1]-psi[2])/(6*grid.dr**2)
    br[0] = 0.
    return br, bz


def azimuthal_current(grid, psi):
    """J_theta = -Delta* psi/(mu0 r), Delta* = d2/dr2 - (1/r) d/dr + d2/dz2 (second differences)."""
    r = grid.r[:, None]
    d2r = np.zeros_like(psi); d2z = np.zeros_like(psi); d1r = np.zeros_like(psi)
    d2r[1:-1] = (psi[2:]-2*psi[1:-1]+psi[:-2])/grid.dr**2
    d1r[1:-1] = (psi[2:]-psi[:-2])/(2*grid.dr)
    d2z[:, 1:-1] = (psi[:, 2:]-2*psi[:, 1:-1]+psi[:, :-2])/grid.dz**2
    with np.errstate(divide="ignore", invalid="ignore"):
        j = -(d2r-d1r/r+d2z)/(MU0*r)
    j[0] = 0.
    return j


def run_display(run, out, label=None, description=""):
    """Write OUT/<run>/frame_NNNN.json (one per saved state) and OUT/<run>/manifest.json (the viewer DATA for this
    run, frames by url); return the run summary for the overview."""
    run = Path(run)
    setup, rf, zf, frames, history = load(run)
    if not setup.get("circuit") or not history[0].get("circuit"):
        raise ValueError(f"{run}: needs a coupled-drive run (drive='coupled') with saved section currents")
    c = case_from_setup(setup)
    i0, r_in, fluid, grid, mapping = geometry(c)
    z0 = mapping.z_offset
    coil = coil_sections(c)
    rr, zz = np.meshgrid(grid.r, grid.z, indexing="ij")
    unit = np.array([f[0] for f in coil.field(rr, zz)])                 # psi per ampere per section
    wall_bz = np.array([float(f[2]) for f in coil.field(c.wall_radius, 0.)])   # vacuum Bz per ampere at the midplane wall
    _, _, _, rows = summarize(run)
    extent = 1.05*max(c.half_length, max(abs(s.z)+s.width/2 for s in coil.straps))
    fluid_rows, fluid_cols = slice(i0, i0+c.nr+1), slice(z0, z0+c.nz+1)
    display_r = grid.r <= c.wall_radius+1.01*fluid.dr                   # keep line colours off the coil filaments
    names = ["Main coil"]*c.sections+[a.name for a in c.aux_coils]
    data_frames, peak_b, peak_j, worst = [], 0., 0., 0.
    (Path(out)/run.name).mkdir(parents=True, exist_ok=True)
    for (k, psi_gas, rho, p), row in zip(frames, rows):
        h = history[k]
        vac = np.asarray(h["circuit"]["vacuum_current_a"]); img = np.asarray(h["circuit"]["image_current_a"])
        source = np.tensordot(vac, unit, axes=1)
        induced = np.zeros(grid.shape)
        induced[fluid_rows, fluid_cols] = psi_gas-source[fluid_rows, fluid_cols]
        psi = mapping.harmonic_extension(induced)+source
        worst = max(worst, float(abs(psi[fluid_rows, fluid_cols]-psi_gas).max()/max(abs(psi_gas).max(), 1e-30)))
        br, bz = fields(grid, psi)
        bmag = np.hypot(br, bz)
        b_at = RegularGridInterpolator((grid.r, grid.z), bmag, bounds_error=False, fill_value=0.)
        gas_b = float(bmag[fluid_rows, fluid_cols].max())
        # Field lines: open flux in equal steps, closed surfaces inside the separatrix.
        region = row if row.get("closed") is not None else {}
        psi_sep = region.get("separatrix_flux_wb", 0.)/(2*np.pi) if region else None
        view = psi[display_r]
        top = float(view.max())
        levels = list(np.linspace(0, top, 13)[1:-1])
        sep_paths = []
        gen = contourpy.contour_generator(x=zz[display_r], y=rr[display_r], z=view)
        if region:
            o = float(psi[fluid_rows, fluid_cols].min())
            levels += list(np.linspace(o, psi_sep, 6)[1:-1])
            lo, hi = zf[max(region["lo"]-1, 0)], zf[min(region["hi"]+1, len(zf)-1)]
            for path in gen.lines(psi_sep-1e-6*abs(psi_sep-o)):
                if path[:, 0].min() >= lo-1e-9 and path[:, 0].max() <= hi+1e-9:
                    sep_paths.extend(rotate(decimate(path), extent))
            if region.get("separatrix_mid_m"):
                sep_paths.append(ring(region["separatrix_mid_m"], 0., extent))
        lines = []
        for level in levels:
            for path in gen.lines(level):
                if len(path) < 3:
                    continue
                path = decimate(path)
                b = [float(f"{x:.3g}") for x in b_at(path[:, ::-1])]
                lines.extend(dict(points=pts, field=b) for pts in rotate(path, extent))
        # Current volume texture: rows r = 0..r_w (support only in the gas), columns over the chamber.
        j = azimuthal_current(grid, psi)[:i0+c.nr+1, fluid_cols]
        support = np.zeros_like(j); support[i0:] = 1.
        j *= support
        texture = np.zeros(j.shape+(2,), dtype="<f4")          # (J_theta, support); the viewer expands to RGBA
        texture[..., 0], texture[..., 1] = j, support
        jmax = float(abs(j).max())
        cloud = []
        sub_r, sub_z = np.arange(i0, i0+c.nr+1, 3), np.arange(0, c.nz+1, 6)   # sparse fallback without WebGL2
        for a in np.arange(6)*np.pi/3:
            for i in sub_r:
                for jj in sub_z:
                    value = float(j[i, jj])
                    if abs(value) > .03*jmax:
                        r = grid.r[i]
                        v = float(f"{value:.3g}")
                        cloud.append([round(zf[jj]/extent, 4), round(r*np.cos(a)/extent, 4), round(r*np.sin(a)/extent, 4), abs(v), 0, v, 0])
        branch = vac+img
        groups = [dict(name=f"{names[s]} {branch[s]/1e3:.1f} kA", center_t=float(branch[s]*wall_bz[s])) for s in range(len(branch))]
        t = float(h["time_s"])
        closed = bool(region.get("closed")) if region else False
        summary = ("No closed flux." if not region else
                   f"{'Closed' if closed else 'Reversed, open at the ends'}: length {100*region['length_m']:.0f} cm · "
                   f"trapped flux {1e3*region['trapped_flux_wb']:.2f} mWb · "
                   + (f"⟨n⟩ {region['density_m3']:.2g} m⁻³ · Te+Ti {region['total_temperature_ev']:.0f} eV" if region.get("density_m3") else ""))
        frame = dict(step=int(h.get("steps", k)), time_us=t*1e6, lines=lines, current=cloud, outline=[], separatrix=sep_paths,
                     bmax=gas_b, jmax=[jmax, 0., jmax, 0.],
                     temperature=[region["total_temperature_ev"]/2]*2 if region.get("total_temperature_ev") else [None, None],
                     drive_center_t=float(sum(g["center_t"] for g in groups)), drive_groups=groups, summary=summary,
                     volume=dict(shape=[texture.shape[1], texture.shape[0]], bounds=[0., c.wall_radius/extent, zf[0]/extent, zf[-1]/extent],
                                 channels=2, slots=[1, 3],    # RGBA = (Jr, Jtheta, Jz, support)
                                 b64=base64.b64encode(texture.tobytes()).decode()))
        if region.get("separatrix_mid_m"):
            frame.update(midplane_radius_m=region["separatrix_mid_m"], maximum_radius_m=region.get("max_separatrix_m"))
        heavy = ("lines", "current", "outline", "separatrix", "volume")      # loaded on demand by the viewer
        name = f"{run.name}/frame_{k:04d}.json"
        (Path(out)/name).write_text(json.dumps({key: frame[key] for key in heavy}, separators=(",", ":")))
        data_frames.append(dict({key: v for key, v in frame.items() if key not in heavy}, url=name))
        peak_b, peak_j = max(peak_b, gas_b), max(peak_j, jmax)
    chamber = [ring(c.wall_radius, z, extent) for z in (-c.half_length, 0., c.half_length)]
    chamber += [(np.array([[-c.half_length, c.wall_radius*np.cos(a), c.wall_radius*np.sin(a)],
                           [c.half_length, c.wall_radius*np.cos(a), c.wall_radius*np.sin(a)]])/extent).round(5).tolist()
                for a in np.arange(8)*np.pi/4]
    if r_in > 0:
        chamber += [ring(r_in, z, extent) for z in (-c.half_length, 0., c.half_length)]
        chamber += [(np.array([[-c.half_length, r_in*np.cos(a), r_in*np.sin(a)], [c.half_length, r_in*np.cos(a), r_in*np.sin(a)]])/extent).round(5).tolist()
                    for a in np.arange(4)*np.pi/2]
    coil_groups = [dict(label=names[s], color=GROUP_COLORS[0] if s < c.sections else GROUP_COLORS[1+(s-c.sections) % (len(GROUP_COLORS)-1)],
                        paths=[ring(st.radius, st.z, extent, 61) for st in coil.straps if st.section == s])
                   for s in range(coil.sections)]
    d = decay(rows)
    title = label or run.name
    gas = f"{c.gas} at {c.fill_pa/0.133322:.1f} mTorr"
    tube = f"{2*c.wall_radius*100/2.54:.0f}\" tube" + (f" with a {2*r_in*100/2.54:.1f}\" inner tube" if r_in > 0 else "")
    scope = (f"Source: {run.name}. {tube}, {gas}, {c.nr}×{c.nz} cells, coupled capacitor banks. Main bank "
             f"{setup['circuit']['main_voltage_v']/1e3:.1f} kV ({setup['circuit']['main_energy_j']/1e3:.1f} kJ). "
             + "".join(f"Extra bank {b['name']}: {b['voltage_v']/1e3:g} kV, {b['energy_j']:.0f} J, fired at {b['fire_time_s']*1e6:g} µs. "
                       for b in setup["circuit"].get("extra_banks", []))
             + f"{len(data_frames)} saved states, every {np.median(np.diff([f['time_us'] for f in data_frames])):.2g} µs. Flux outside the gas is the run's own vacuum solve from the saved "
             f"gas-node flux and section currents (gas-node flux reproduced to {worst:.0e}). Field lines are flux contours "
             "(equal flux between open lines, closed surfaces inside the gold separatrix), coloured by |B|, with the colour "
             "scale set by the field in the gas. Current is J_θ = −Δ*ψ/(μ0 r); J_r and J_z are zero because the model has no "
             "toroidal field. 2D axisymmetric single-fluid MHD: it does not model the n = 2 rotational mode, so FRC lifetime "
             "beyond the simulated window is not predicted."
             + (f" Trapped-flux e-fold time fitted from 8 µs: {d['flux_efold_s']*1e6:.0f} µs." if d.get("flux_efold_s") else ""))
    data = dict(frames=data_frames, detailed=True, bmax=peak_b, jmax=[peak_j, 0., peak_j, 0.],
                initial_frame=len(data_frames)-1, playback_fps=20, initial_yaw=-.35, initial_pitch=.28, initial_zoom=3.,
                title=title, description=f"{description or 'Saved FRC formation run'}. Gold: separatrix. Copper and coloured rings: coils.",
                chamber=chamber, coils=[], coil_groups=coil_groups, physical_scale_m=extent,
                status_note=f"{title} · saved formation states", scope=scope, temperature_digits=3,
                drive_label="Vacuum Bz of the coil currents at the tube wall, midplane", motion_label="Separatrix radius — midplane",
                solver_steps=int(history[-1].get("steps", 0)))
    particles = particle_trails(run, out, extent)
    if particles:
        data["particles"] = particles
    (Path(out)/run.name/"manifest.json").write_text(json.dumps(data, separators=(",", ":")))
    return dict(run=run.name, label=title, description=description, manifest=f"{run.name}/manifest.json", frames=len(data_frames), flux_reproduction=worst,
                particles=bool(particles),
                inner_radius_m=r_in, gas=gas, banks=scope,
                rows=[dict({k: r.get(k) for k in ("time_s", "closed", "separatrix_mid_m", "length_m", "trapped_flux_wb", "density_m3",
                                                  "total_temperature_ev")}, formed=formed(r, 2*c.half_length)) for r in rows],
                flux_efold_s=d.get("flux_efold_s"), fit_span_s=(d.get("last_closed_s") or 0.)-8e-6)


def viewer_page(out, summaries):
    """OUT/viewer.html: one viewer for all runs. It embeds the first run's manifest and a run selector; switching
    fetches the other manifests and keeps the camera and the simulated time (viewer.js switchRun; ?run=NAME opens one)."""
    first = json.loads((Path(out)/summaries[0]["manifest"]).read_text())
    data = dict(first, runs=[dict(name=s["run"], label=s["label"], manifest=s["manifest"]) for s in summaries], run=summaries[0]["run"])
    html = (ASSETS/"volume.html").read_text()
    reps = [("<title>FRC sandbox · Plasma volume</title>", "<title>FRC formation runs</title>"),
            ("<h1>Inside the plasma</h1>", '<h1 id="run-title"></h1>'),
            ("<p>Current density in volume. Magnetic structure in motion. Two views of the same saved state.</p>", '<p id="run-description"></p>'),
            ('<div class="toolbar">', '<div class="toolbar"><label for="run">Run</label><select id="run" aria-label="Run"></select>'),
            ("Axisymmetric · 3D rendering", "2D axisymmetric runs · 3D rendering"),
            ('<p id="motion"></p>', '<p id="motion"></p><p id="frame-summary"></p>'),
            ('<label><input type="checkbox" id="geometry"> Tube & coils</label>',
             '<label><input type="checkbox" id="geometry"> Tube & coils</label><label><input type="checkbox" id="separatrix" checked> Separatrix</label>'),
            ("Mean electron / ion temperature", "Te / Ti inside separatrix (Te = Ti)"),
            (">96<", ">×1<"), (">192<", ">×2<"), (">384<", ">×4<"),
            # No toroidal field in the model, so Ampère gives J_r = J_z = 0 exactly: shown, but not selectable.
            ('<option value="1">Radial Jr</option>', '<option value="1" disabled>Radial Jr (zero: no toroidal field in this model)</option>'),
            ('<option value="3">Axial Jz</option>', '<option value="3" disabled>Axial Jz (zero: no toroidal field in this model)</option>'),
            ("Exploratory pre-ionized plasma model. More rendering detail does not add simulation resolution or establish FRC formation, stability or confinement.",
             "Assumed complete pre-ionization; 2D axisymmetric MHD. Rendering detail adds no simulation resolution; the n = 2 mode and lifetime are not modelled.")]
    if any(s.get("particles") for s in summaries):
        reps += [('<label><input type="checkbox" id="separatrix" checked> Separatrix</label>',
                  '<label><input type="checkbox" id="separatrix" checked> Separatrix</label><label><input type="checkbox" id="particles"> Particles</label>'
                  '<select id="particle-set" aria-label="Particle orbit set" hidden></select>'),
                 ('<p id="frame-summary"></p>', '<p id="frame-summary"></p><p id="particle-legend"></p>')]
    for a, b in reps:
        if a not in html:
            raise ValueError(f"volume.html no longer contains {a!r}")
        html = html.replace(a, b, 1)
    html = html.replace("/* DATA */", "const DATA="+json.dumps(data, separators=(",", ":"))+";")
    html = html.replace("/* VOLUME */", (ASSETS/"volume.js").read_text()).replace("/* VIEWER */", (ASSETS/"viewer.js").read_text())
    (Path(out)/"viewer.html").write_text(html)


def describe(setup):
    """One-line description of an experiment from its saved case (what differs from a plain FRTP)."""
    c = case_from_setup(setup)
    coil = "uniform" if not c.taper_start else f"graded {100*c.coil_radius:.1f}→{100*c.end_radius:.1f} cm"
    text = [f"{c.gas} {c.fill_pa/0.133322:.1f} mTorr, {coil} coil, bias {c.bias_t:+.2f} T, peak {c.peak_t:.2f} T "
            f"({setup['circuit']['main_voltage_v']/1e3:.1f} kV, {setup['circuit']['main_energy_j']/1e3:.1f} kJ)"]
    i0, r_in = inner_rows(c)
    if r_in > 0:
        text.append(f"inner tube r {100*r_in:.2f} cm")
    for b in c.extra_main_banks:
        text.append(f"re-boost bank {b.capacitance*1e6:g} µF at {b.voltage/1e3:g} kV fired at {b.fire_time*1e6:g} µs")
    for a in c.aux_coils:
        where = f"z ±{a.z_center:g} m" if a.pair else f"z {a.z_center:g} m"
        text.append(f"{a.name}: {a.turns}-turn coil r {100*a.radius:g} cm × {100*a.length:g} cm at {where}, "
                    f"{a.bank.capacitance*1e6:g} µF at {a.bank.voltage/1e3:g} kV fired at {a.bank.fire_time*1e6:g} µs")
    return "; ".join(text)


OVERVIEW = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Formation runs</title>
<style>
:root{color-scheme:dark;--surface-0:#080b0f;--surface-1:#10161e;--border:#26303b;--grid:#1f2833;--text-primary:#e8edf2;--text-secondary:#98a6b5;--text-muted:#6f7f8d;--accent:#7ce5d2;
--series-1:#3987e5;--series-2:#d95926;--series-3:#199e70;--series-4:#c98500;--series-5:#d55181;--series-6:#008300}
*{box-sizing:border-box}body{margin:0;background:var(--surface-0);color:var(--text-primary);font:14px system-ui,sans-serif}
main{max-width:1280px;margin:auto;padding:26px 16px 40px}h1{font-weight:500;font-size:30px;letter-spacing:-.03em;margin:8px 0}
p{line-height:1.6;color:var(--text-secondary);max-width:900px}.eyebrow{color:var(--accent);font-size:10px;letter-spacing:.19em;text-transform:uppercase}
.legend{display:flex;flex-wrap:wrap;gap:6px 14px;margin:18px 0 6px}.legend button{display:flex;align-items:center;gap:7px;background:none;border:1px solid var(--border);border-radius:6px;color:var(--text-primary);padding:6px 9px;font:inherit;font-size:12px;cursor:pointer}
.legend button[aria-pressed=false]{color:var(--text-muted);opacity:.6}.swatch{width:18px;height:0;border-top:2px solid}.swatch.dashed{border-top-style:dashed}
.charts{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.chart{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;padding:12px 12px 6px;position:relative}
.chart h2{font-size:14px;font-weight:500;margin:0 0 2px}.chart small{color:var(--text-muted);font-size:11px}svg{display:block;width:100%;height:auto;overflow:visible}
.axis text{fill:var(--text-muted);font-size:10px}.axis line,.grid line{stroke:var(--grid)}.tip{position:absolute;pointer-events:none;background:#0b1016;border:1px solid var(--border);border-radius:6px;padding:7px 9px;font-size:11px;color:var(--text-primary);display:none;white-space:nowrap;z-index:3}
.tip div{display:flex;align-items:center;gap:6px;color:var(--text-secondary)}.tip b{color:var(--text-primary);font-weight:500}
table{width:100%;border-collapse:collapse;margin-top:22px;font-size:12px}th,td{text-align:left;padding:8px 9px;border-bottom:1px solid var(--border);vertical-align:top}th{color:var(--text-muted);font-weight:500}
td{color:var(--text-secondary)}td.name{color:var(--text-primary);white-space:nowrap}td.num{font-variant-numeric:tabular-nums;white-space:nowrap}a{color:var(--accent)}
.note{font-size:12px;color:var(--text-muted)}@media(max-width:900px){.charts{grid-template-columns:1fr}table{display:block;overflow-x:auto}}
</style>
<main><div class="eyebrow">FRC experiment playground · 12-inch device</div><h1>Formation runs compared</h1>
<p><a href=\"viewer.html\">Open the 3D viewer</a> (choose a run in it). Each line is one simulated experiment. Values are drawn only while a formed FRC exists: closed flux that touches no wall and is shorter than 80% of the tube (the design-scan criterion); a gap means it is not. Solid lines: no inner tube; dashed: 3" inner tube. Hover a chart for exact values; open a run's 3D view from the table.</p>
<div class="legend" id="legend"></div>
<div class="charts" id="charts"></div>
<table id="table"></table>
<p class="note" id="scope"></p></main>
<script>
const DATA=/* DATA */;
const charts=[{key:'trapped_flux_wb',title:'Trapped (closed) flux',unit:'mWb',scale:1e3,digits:2},{key:'length_m',title:'Separatrix length',unit:'cm',scale:100,digits:0},{key:'separatrix_mid_m',title:'Separatrix radius, midplane',unit:'cm',scale:100,digits:1}];
const visible=DATA.runs.map(()=>true),tmax=Math.max(...DATA.runs.flatMap(r=>r.rows.map(x=>x.time_s*1e6)));
const el=(tag,attrs={},text)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const k in attrs)e.setAttribute(k,attrs[k]);if(text!=null)e.textContent=text;return e};
const value=(row,c)=>row.formed&&row[c.key]!=null?row[c.key]*c.scale:null;
function nearest(r,t){let best=null,gap=Infinity;for(const q of r.rows){const g=Math.abs(q.time_s*1e6-t);if(g<gap){gap=g;best=q}}return gap<=.55?best:null}
function nice(max){const p=Math.pow(10,Math.floor(Math.log10(max))),m=max/p;return [1,1.2,1.6,2,2.4,3,4,5,6,8,10].find(k=>k>=m-1e-9)*p}
function legend(){const L=document.getElementById('legend');L.replaceChildren();DATA.runs.forEach((r,i)=>{const b=document.createElement('button');b.setAttribute('aria-pressed',visible[i]);
 const s=document.createElement('span');s.className='swatch'+(r.dashed?' dashed':'');s.style.borderTopColor=`var(--series-${i+1})`;b.append(s,document.createTextNode(r.label));
 b.onclick=()=>{visible[i]=!visible[i];draw()};L.append(b)})}
function draw(){legend();const C=document.getElementById('charts');C.replaceChildren();
 charts.forEach(c=>{const box=document.createElement('section');box.className='chart';const h=document.createElement('h2');h.textContent=c.title;const u=document.createElement('small');u.textContent=c.unit+' vs time (µs)';box.append(h,u);
  const W=380,H=230,m={l:34,r:8,t:10,b:22},x=t=>m.l+(W-m.l-m.r)*t/tmax;const vals=DATA.runs.flatMap((r,i)=>visible[i]?r.rows.map(row=>value(row,c)).filter(v=>v!=null):[]);
  const ymax=nice(Math.max(1e-9,...vals)*1.05),y=v=>H-m.b-(H-m.t-m.b)*v/ymax;const svg=el('svg',{viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':c.title});
  const g=el('g',{class:'grid'}),a=el('g',{class:'axis'});for(let k=0;k<=4;k++){const v=ymax*k/4;g.append(el('line',{x1:m.l,x2:W-m.r,y1:y(v),y2:y(v)}));a.append(el('text',{x:m.l-5,y:y(v)+3,'text-anchor':'end'},+(v).toPrecision(3)))}
  const step=tmax>30?10:5;for(let t=0;t<=tmax+1e-9;t+=step)a.append(el('text',{x:x(t),y:H-6,'text-anchor':'middle'},t));svg.append(g,a);
  DATA.runs.forEach((r,i)=>{if(!visible[i])return;let d='',pen=false;r.rows.forEach(row=>{const v=value(row,c);if(v==null){pen=false;return}d+=(pen?'L':'M')+x(row.time_s*1e6).toFixed(1)+' '+y(v).toFixed(1);pen=true});
   svg.append(el('path',{d,fill:'none',stroke:`var(--series-${i+1})`,'stroke-width':2,'stroke-linejoin':'round','stroke-linecap':'round','stroke-dasharray':r.dashed?'5 3':'none'}))});
  const cross=el('line',{y1:m.t,y2:H-m.b,stroke:'#5b6b78','stroke-width':1,visibility:'hidden'});svg.append(cross);
  const hit=el('rect',{x:m.l,y:m.t,width:W-m.l-m.r,height:H-m.t-m.b,fill:'transparent'});svg.append(hit);const tip=document.createElement('div');tip.className='tip';box.append(svg,tip);
  hit.addEventListener('pointermove',e=>{const rect=svg.getBoundingClientRect(),sx=(e.clientX-rect.left)*W/rect.width,t=Math.round(10*Math.max(0,Math.min(tmax,(sx-m.l)/(W-m.l-m.r)*tmax)))/10;
   cross.setAttribute('x1',x(t));cross.setAttribute('x2',x(t));cross.setAttribute('visibility','visible');tip.replaceChildren();const head=document.createElement('b');head.textContent=t.toFixed(1)+' µs';tip.append(head);
   DATA.runs.forEach((r,i)=>{if(!visible[i])return;const row=nearest(r,t);const v=row?value(row,c):null;const line=document.createElement('div');const s=document.createElement('span');s.className='swatch'+(r.dashed?' dashed':'');s.style.borderTopColor=`var(--series-${i+1})`;
    line.append(s,document.createTextNode(r.label+': '+(v==null?(row?'open':'—'):v.toFixed(c.digits)+' '+c.unit)));tip.append(line)});
   tip.style.display='block';const px=(e.clientX-box.getBoundingClientRect().left);tip.style.left=Math.min(px+12,box.clientWidth-tip.offsetWidth-4)+'px';tip.style.top='40px'});
  hit.addEventListener('pointerleave',()=>{tip.style.display='none';cross.setAttribute('visibility','hidden')});C.append(box)})}
function table(){const T=document.getElementById('table');const at=(r,t,k,s,d)=>{const row=nearest(r,t);return row&&row.formed&&row[k]!=null?(row[k]*s).toFixed(d):row?'open':'—'};
 const times=[10,20,30,40].filter(t=>t<=tmax+.5);T.innerHTML='';const head=T.insertRow();['Experiment','Setup','Formed from (held 3 µs)']
  .concat(times.map(t=>'φ at '+t+' µs (mWb)'),['l_s at 20 µs (cm)','Flux e-fold','3D view']).forEach(h=>{const th=document.createElement('th');th.textContent=h;head.append(th)});
 DATA.runs.forEach((r,i)=>{const tr=T.insertRow();const name=tr.insertCell();name.className='name';const s=document.createElement('span');s.className='swatch'+(r.dashed?' dashed':'');s.style.cssText=`display:inline-block;vertical-align:middle;margin-right:7px;border-top-color:var(--series-${i+1})`;name.append(s,document.createTextNode(r.label));
  tr.insertCell().textContent=r.description;const first=r.rows.find(q=>r.rows.every(o=>o.time_s<q.time_s-1e-9||o.time_s>q.time_s+3e-6||o.formed));const c0=tr.insertCell();c0.className='num';c0.textContent=first?(+(first.time_s*1e6).toFixed(1))+' µs':'never';
  times.forEach(t=>{const c=tr.insertCell();c.className='num';c.textContent=at(r,t,'trapped_flux_wb',1e3,2)});const ls=tr.insertCell();ls.className='num';ls.textContent=at(r,20,'length_m',100,0);
  const ef=tr.insertCell();ef.className='num';ef.textContent=r.flux_efold_s?Math.round(r.flux_efold_s*1e6)+' µs':'—';const link=tr.insertCell();const a=document.createElement('a');a.href=r.page;a.textContent='Open';link.append(a)})}
document.getElementById('scope').textContent=DATA.scope;draw();table();
</script></html>"""


def overview(out, summaries):
    build = int(time.time())             # version tag on links, so a rebuilt page is not served from the browser cache
    data = dict(runs=[dict(label=s["label"], description=s["description"], page=f"viewer.html?run={s['run']}&v={build}", dashed=s["inner_radius_m"] > 0,
                           rows=s["rows"], flux_efold_s=s["flux_efold_s"] if s["fit_span_s"] >= 15e-6 else None) for s in summaries],
                scope="2D axisymmetric single-fluid MHD with coupled capacitor banks; assumed complete pre-ionization. "
                      "The n = 2 rotational mode is not modelled, so lifetime beyond the simulated window is not predicted. "
                      "Flux e-fold: least-squares fit of ln φ over closed frames from 8 µs, shown when the fit spans at least 15 µs. "
                      "Generated by examples/formation_display.py.")
    if len(summaries) > 6:
        raise ValueError("The overview's validated categorical palette has six slots")
    (Path(out)/"index.html").write_text(OVERVIEW.replace("/* DATA */", json.dumps(data, default=float)))


def main(out, runs):
    """runs: RUN or RUN::Label. Writes OUT/viewer.html (all runs, run selector), OUT/<run>/ (manifest and frames)
    and OUT/index.html (overview)."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    summaries = []
    for spec in runs:
        run, _, label = spec.partition("::")
        setup = json.loads((Path(run)/"setup.json").read_text())
        summaries.append(run_display(run, out, label or None, describe(setup)))
        print(f"display {summaries[-1]['run']}: {summaries[-1]['frames']} frames, gas flux reproduced to {summaries[-1]['flux_reproduction']:.1e}", flush=True)
    (out/"runs.json").write_text(json.dumps(summaries, indent=1, default=float))
    for stale in out.glob("*.html"):           # per-run pages of earlier builds
        stale.unlink()
    viewer_page(out, summaries)
    overview(out, summaries)
    return summaries


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
