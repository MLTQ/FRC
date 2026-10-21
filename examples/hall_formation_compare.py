"""Compare a Hall (2.5D) formation run with a poloidal-only reference of the same experiment.

`python -m examples.hall_formation_compare HALL_RUN REFERENCE_RUN [OUT.json]`

Both runs come from examples/nucte_formation.py. Frames are matched by saved time. Reported:
- end reconnection: the first time from which the FRC stays formed (closed, shorter than 80% of the
  tube) for 3 us, as in the design scan, and the first times the closed length falls below 80%, 60%
  and 40% of the tube (the separatrix ends pulling in from the tube ends);
- contraction: the time and value of the minimum separatrix length after closure;
- r_s (midplane), l_s, trapped flux, <n> and Te+Ti at fixed times;
- the Hall run's peak B_theta, v_theta, J_r, J_z (history rows) and the B_theta energy;
- steps and wall time (wall time depends on machine load; steps do not);
- for the Hall run, peak |v_theta| and |B_theta| and the v_theta kinetic-energy share by density band
  (n/fill below 0.05: near the 1% gap floor; 0.05-0.3; 0.3-1; above 1: the FRC plasma), every 5 us.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.constants import k as KB
from examples.nucte_compare import summarize, load
from examples.nucte_formation import ION_MASS
from examples.design_scan_summary import formed

TIMES = (2e-6, 4e-6, 6e-6, 8e-6, 10e-6, 12e-6, 15e-6, 20e-6)


def timeline(run):
    setup, rf, zf, rows = summarize(run)
    history = json.loads((Path(run)/"history.json").read_text())
    tube = float(zf[-1]-zf[0])
    t = np.array([r["time_s"] for r in rows])
    ok = np.array([formed(r, tube) for r in rows])
    closed = np.array([bool(r.get("closed")) for r in rows])
    first_closed = float(t[closed][0]) if closed.any() else None
    held = next((float(t[k]) for k in range(len(t)) if ok[k] and ok[(t >= t[k]) & (t <= t[k]+3e-6)].all()), None)
    lengths = np.array([r.get("length_m", np.nan) if r.get("closed") else np.nan for r in rows])
    below = {f"length_below_{int(100*f)}pct_s": next((float(t[k]) for k in range(len(t)) if lengths[k] < f*tube), None)
             for f in (.8, .6, .4)}
    after = closed & (t >= (first_closed or 0.))
    k_min = int(np.nanargmin(np.where(after, lengths, np.nan))) if after.any() else None
    at = {}
    for time in TIMES:
        k = int(np.argmin(abs(t-time)))
        if abs(t[k]-time) > 5e-8:
            continue
        r = rows[k]
        at[f"{time*1e6:g}us"] = dict(closed=bool(r.get("closed")), formed=bool(ok[k]),
                                     separatrix_mid_cm=None if r.get("separatrix_mid_m") is None else 100*r["separatrix_mid_m"],
                                     length_cm=None if not r.get("closed") else 100*r["length_m"],
                                     trapped_flux_mwb=None if not r.get("closed") else 1e3*r["trapped_flux_wb"],
                                     density_m3=r.get("density_m3"), te_plus_ti_ev=r.get("total_temperature_ev"))
    last = history[-1]
    peaks = {key: max((h.get(key) or 0.) for h in history)
             for key in ("max_btheta_t", "max_vtheta_m_s", "max_jz_a_m2", "max_jr_a_m2", "toroidal_magnetic_energy_j",
                         "toroidal_kinetic_energy_j")}
    return dict(run=str(run), end_s=float(t[-1]), first_closed_s=first_closed, formed_held_from_s=held,
                minimum_length_s=None if k_min is None else float(t[k_min]),
                minimum_length_cm=None if k_min is None else 100*float(lengths[k_min]), **below, at=at,
                steps=last.get("steps"), wall_s=last.get("wall_s"),
                response_energy_balance_relative=last.get("response_energy_balance_relative"),
                hall_work_j=last.get("hall_work_j"), hall_limited_steps=last.get("hall_limited_steps"),
                peaks=peaks if any(peaks.values()) else None)


def toroidal_by_density(run, every_s=5e-6):
    """Where the toroidal velocity and field live, by density band, from toroidal.npz and the snapshots."""
    path = Path(run)/"toroidal.npz"
    if not path.exists():
        return None
    setup, rf, zf, frames, history = load(run)
    case = setup["case"]
    fill = 2*case["fill_pa"]/(KB*case["gas_k"])
    data, nr, nz = np.load(path), len(rf)-1, len(zf)-1
    radius = np.sqrt(.5*(rf[:-1]**2+rf[1:]**2))[:, None]
    bands = ((0., .05), (.05, .3), (.3, 1.), (1., np.inf))
    out = {}
    for k, psi, rho, p in frames:
        time = history[k]["time_s"]
        if k == 0 or abs(time/every_s-round(time/every_s)) > 1e-3 or f"t{k:04d}" not in data.files:
            continue
        v = data[f"t{k:04d}"]
        vt, bt = radius*v[:nr*nz].reshape(nr, nz), radius*v[nr*nz:].reshape(nr, nz)
        n = rho/ION_MASS[case["gas"]]/fill
        kinetic = .5*rho*vt**2*radius
        row = {}
        for lo, hi in bands:
            m = (n >= lo) & (n < hi)
            row[f"n/fill {lo:g}-{hi:g}"] = dict(cells=int(m.sum()), vtheta_energy_share=float(kinetic[m].sum()/kinetic.sum()),
                                               max_vtheta_km_s=float(abs(vt[m]).max()/1e3) if m.any() else 0.,
                                               max_btheta_mt=float(abs(bt[m]).max()*1e3) if m.any() else 0.)
        out[f"{time*1e6:g}us"] = row
    return out


def steps_at(run, time):
    history = json.loads((Path(run)/"history.json").read_text())
    row = min(history, key=lambda h: abs(h["time_s"]-time))
    return row["steps"], row.get("wall_s")


def main(hall, reference, out=None):
    a, b = timeline(hall), timeline(reference)
    end = min(a["end_s"], b["end_s"])
    steps = dict(hall=steps_at(hall, end), reference=steps_at(reference, end))
    result = dict(hall=a, reference=b, common_end_s=end, steps_and_wall_at_common_end=steps,
                  step_ratio=steps["hall"][0]/steps["reference"][0], hall_toroidal_by_density=toroidal_by_density(hall))
    text = json.dumps(result, indent=1, default=float)
    if out:
        Path(out).write_text(text)
    print(text)
    return result


if __name__ == "__main__":
    main(*sys.argv[1:4])
