"""Summarize formation runs of a design scan: bank numbers, when a closed FRC forms, and its state."""
import json
import sys
from pathlib import Path
import numpy as np
from examples.nucte_compare import summarize, decay, coil_profile


def formed(row, chamber_length):
    """Closed FRC: separatrix off the wall and its ends clear of the chamber ends (< 80% of its length)."""
    return bool(row.get("closed")) and row.get("length_m", np.inf) < .8*chamber_length


def run_summary(run, times=(10e-6, 15e-6, 20e-6)):
    setup, rf, zf, rows = summarize(run)
    case = setup["case"]
    circuit = setup.get("circuit") or {}
    chamber = float(zf[-1]-zf[0])
    ok = [formed(r, chamber) for r in rows]
    t = [r["time_s"] for r in rows]
    # first time from which the FRC stays closed for 3 us (or to the end of the run)
    first = next((t[k] for k in range(len(rows)) if all(ok[j] for j in range(k, len(rows)) if t[j] <= t[k]+3e-6)), None)
    out = dict(run=Path(run).name, gas=case.get("gas", "D2"), fill_mtorr=case["fill_pa"]/0.133322, bias_t=case["bias_t"],
               peak_t=case["peak_t"], coil=f"{case['coil_radius']:.3f}->{case['end_radius']:.3f}" if case.get("taper_start") else f"{case['coil_radius']:.3f}",
               charge_kv=circuit.get("main_voltage_v", np.nan)/1e3, bank_kj=circuit.get("main_energy_j", np.nan)/1e3,
               formed_s=first, last_time_s=rows[-1]["time_s"], decay=decay(rows), coil_shape=coil_profile(case))
    for t in times:
        r = min(rows, key=lambda x: abs(x["time_s"]-t))
        if abs(r["time_s"]-t) > .6e-6:
            continue
        out[f"{t*1e6:g}us"] = dict(formed=formed(r, chamber), rs_cm=None if r.get("separatrix_mid_m") is None else 100*r["separatrix_mid_m"],
                                   ls_cm=None if r.get("length_m") is None else 100*r["length_m"], n_m3=r.get("density_m3"),
                                   t_ev=r.get("total_temperature_ev"), phi_mwb=None if r.get("trapped_flux_wb") is None else 1e3*r["trapped_flux_wb"])
    return out


def fmt(v, spec):
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, spec)


def main(runs):
    rows = [run_summary(r) for r in runs]
    print("| Run | Gas, fill | Bias, peak | Coil r (m) | R_M | Charge, bank | Closed FRC from | At 10 µs: r_s, l_s, n, T, φ | At 15–20 µs | Flux e-fold |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for s in rows:
        def state(key):
            d = s.get(key)
            if not d:
                return "—"
            tag = "" if d["formed"] else " (open)"
            return f"{fmt(d['rs_cm'], '.1f')} cm, {fmt(d['ls_cm'], '.0f')} cm, {fmt(d['n_m3'], '.2g')}, {fmt(d['t_ev'], '.0f')} eV, {fmt(d['phi_mwb'], '.2f')} mWb{tag}"
        late = s.get("20us") or s.get("15us")
        late_key = "20us" if s.get("20us") else "15us"
        print(f"| {s['run']} | {s['gas']}, {s['fill_mtorr']:.1f} mTorr | {s['bias_t']:+.2f} T, {s['peak_t']:.2f} T | {s['coil']} | "
              f"{s['coil_shape']['mirror_ratio']:.2f} | {fmt(s['charge_kv'], '.1f')} kV, {fmt(s['bank_kj'], '.1f')} kJ | "
              f"{fmt(None if s['formed_s'] is None else s['formed_s']*1e6, '.0f')} µs | {state('10us')} | {state(late_key)} | "
              f"{fmt(None if s['decay']['flux_efold_s'] is None else s['decay']['flux_efold_s']*1e6, '.0f')} µs |")
    Path("outputs/design_scan").mkdir(parents=True, exist_ok=True)
    Path("outputs/design_scan/summary.json").write_text(json.dumps(rows, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1:])


def signatures(run):
    """Per-microsecond signals the lab diagnostics would record (midplane unless stated):
    wall Bz (B-dot probe just inside the tube), excluded-flux radius r_dphi from a flux loop at the wall,
    side-on interferometer integral of n dl through the axis, separatrix length (camera), and the bank's
    coil current and module capacitor voltage (Rogowski, voltage divider)."""
    import numpy as np
    from examples.nucte_compare import load
    setup, rf, zf, frames, history = load(run)
    mass = {"D2": 2.01410, "H2": 1.00784}[setup["case"].get("gas", "D2")]*1.66053906660e-27
    j = len(zf)//2
    out = []
    for k, psi, rho, p in frames:
        h = history[k] if k < len(history) else {}
        b_wall = (psi[-1, j]-psi[-2, j])/(rf[-1]-rf[-2])*2/(rf[-1]+rf[-2])
        r_dphi = np.sqrt(max(rf[-1]**2-2*psi[-1, j]/b_wall, 0.)) if b_wall != 0 else np.nan
        chord = float(2*np.sum(.5*(rho[:, j-1]+rho[:, j])/mass*np.diff(rf)))
        c = h.get("circuit") or {}
        coil = sum(c.get("vacuum_current_a", [np.nan]))+sum(c.get("image_current_a", [0.]))
        caps = [v for v, m in zip(c.get("capacitor_voltage_v", []), range(10**6))]
        out.append(dict(time_us=h.get("time_s", np.nan)*1e6, wall_bz_t=float(b_wall), excluded_flux_radius_cm=100*float(r_dphi),
                        interferometer_ndl_m2=chord, separated_length_m=h.get("separated_length_m"),
                        coil_current_ka=coil/1e3, main_capacitor_kv=(caps[-1]/1e3 if caps else np.nan)))
    return out
