"""First-cut benchtop FRC design window: independent formula checks, scan and sensitivity."""
import csv
import json
import sys
from dataclasses import replace, asdict
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from plasma_bench.electromagnetics import loop_field
from plasma_bench.frc_scaling import (Assumptions, Device, design_point, rigid_rotor_flux,
                                      rigid_rotor_decay, sheet_coil_center_field,
                                      wheeler_inductance, record)

GASES = ("H2", "D2", "N2")
FILLS_MTORR = (3., 5., 10., 20., 40.)
FIELDS_T = (.1, .2, .3, .5, .75)
WINDOW = dict(s_min=1., s_max=4., te_min_ev=10., energy_max_j=5e3, voltage_max_v=25e3)
RANGES = dict(particle_retention=(.2, .6), flux_retention=(.2, .6), separatrix_to_wall=(.4, .6),
              length_to_coil=(.4, .8), shape=(2/3, .9), profile_k=(1., 3.), ti_over_te=(.5, 2.))
KEYS = ("gas", "fill_mtorr", "b_ext_t", "total_t_ev", "te_ev", "n_avg_m3", "s_bar", "rs_over_di",
        "tau_phi_classical_s", "bias_t", "trapped_flux_wb", "snowplow_s", "quarter_period_s",
        "coil_current_a", "turns", "segments", "module_current_a", "loop_voltage_v",
        "bank_energy_j", "bank_voltage_v", "module_capacitance_f", "bank_capacitance_f", "pi_te_min_ev")
DRIVES = ((1, 1), (2, 1), (3, 1), (4, 1), (1, 2), (1, 4), (2, 2), (1, 8))  # (series turns, parallel modules)


def check_rigid_rotor(k, r_s=1., b_e=1.):
    """Quadrature of the tanh profile and a finite-difference null current, independent of closed forms."""
    bz = lambda r: b_e*np.tanh(k*(2*r*r/r_s**2-1))
    null = r_s/np.sqrt(2)
    flux = quad(lambda r: bz(r)*2*np.pi*r, 0, null, epsabs=0, epsrel=1e-12)[0]
    sep = quad(lambda r: bz(r)*2*np.pi*r, 0, r_s, epsabs=1e-12, epsrel=1e-12)[0]
    h = 1e-6*r_s
    dbdr = (bz(null+h)-bz(null-h))/(2*h)              # = mu_0 J_theta at the null
    tau_eta = abs(flux)/(2*np.pi*null*dbdr*r_s**2)     # tau*eta/(mu_0 r_s^2), eta J = E_theta
    return dict(k=k, flux_quadrature=abs(flux)/(np.pi*r_s**2*b_e), flux_closed=float(rigid_rotor_flux(k)),
                separatrix_flux=sep/(np.pi*r_s**2*b_e), decay_numeric=tau_eta,
                decay_closed=float(rigid_rotor_decay(k)))


def check_coil(radius, length, rings=4000):
    """Sum analytic elliptic-integral loops for a uniform sheet; compare centre field and inductance."""
    z_src = (np.arange(rings)+.5)/rings*length-length/2
    _, _, bz0 = loop_field(np.zeros(rings), -z_src, radius)
    b_centre = bz0.sum()/rings
    z_mid = z_src[:-1]+length/rings/2
    psi = np.zeros(z_mid.size)
    for z0 in z_src:
        psi += loop_field(np.full(z_mid.size, radius), z_mid, radius, z0)[0]/rings
    inductance = 2*np.pi*psi.mean()*(rings-1)/rings+2*np.pi*psi[[0, -1]].mean()/rings
    return dict(radius_m=radius, length_m=length, rings=rings,
                centre_field_sum_t_per_a=float(b_centre),
                centre_field_formula_t_per_a=float(sheet_coil_center_field(1., radius, length)),
                inductance_flux_linkage_h=float(inductance),
                inductance_wheeler_h=float(wheeler_inductance(radius, length)))


def in_window(p):
    return (WINDOW["s_min"] <= p["s_bar"] <= WINDOW["s_max"] and p["te_ev"] >= WINDOW["te_min_ev"]
            and p["bank_energy_j"] <= WINDOW["energy_max_j"] and p["bank_voltage_v"] <= WINDOW["voltage_max_v"])


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(out=Path("outputs/frc_design_window"), gas="D2", fill=10., field=.3):
    out.mkdir(parents=True, exist_ok=True)
    device, base = Device(), Assumptions()
    checks = dict(rigid_rotor=[check_rigid_rotor(k) for k in (1., 2., 3.)],
                  coil=[check_coil(device.coil_radius, device.coil_length, n) for n in (1000, 4000)])
    scan = []
    for g in GASES:
        for p in FILLS_MTORR:
            for b in FIELDS_T:
                row = {k: design_point(g, p, b, device, base)[k] for k in KEYS}
                scan.append(dict(row, in_window=in_window(row)))
    sensitivity = [dict(case="nominal", **{k: design_point(gas, fill, field, device, base)[k] for k in KEYS})]
    for name, values in RANGES.items():
        for v in values:
            a = replace(base, **{name: v})
            sensitivity.append(dict(case=f"{name}={v:.3g}", **{k: design_point(gas, fill, field, device, a)[k] for k in KEYS}))
    for v in (15e-9, 60e-9):
        sensitivity.append(dict(case=f"stray_inductance={v:.3g}",
                                **{k: design_point(gas, fill, field, replace(device, stray_inductance=v), base)[k] for k in KEYS}))
    for r_w in (.05, .10):
        d = replace(device, wall_radius=r_w, coil_radius=r_w+.009)
        sensitivity.append(dict(case=f"wall_radius={r_w:.3g}", **{k: design_point(gas, fill, field, d, base)[k] for k in KEYS}))
    drives = [{k: design_point(g, fill, field, replace(device, turns=n, segments=m), base)[k] for k in KEYS}
              for g in ("H2", "D2") for n, m in DRIVES]
    write_csv(out/"drive_options.csv", drives)
    write_csv(out/"scan.csv", scan)
    write_csv(out/"sensitivity.csv", sensitivity)
    (out/"checks.json").write_text(json.dumps(checks, indent=2))
    (out/"summary.json").write_text(json.dumps(dict(record(device, base), window=WINDOW, ranges=RANGES,
                                                    nominal=dict(gas=gas, fill_mtorr=fill, b_ext_t=field)), indent=2))
    json.dump(checks, sys.stdout, indent=1)
    print()
    for row in sensitivity:
        print("{case:>26} T={total_t_ev:6.1f}eV n={n_avg_m3:.2e} s={s_bar:4.2f} tau={tau_phi_classical_s:.2e}s "
              "bias={bias_t:+.4f}T E={bank_energy_j:6.0f}J V={bank_voltage_v:6.0f}V".format(**row))
    for row in drives:
        print("{gas} turns={turns} modules={segments}: {module_current_a:8.0f} A/module {bank_voltage_v:6.0f} V "
              "{module_capacitance_f:.2e} F/module total {bank_energy_j:4.0f} J, loop {loop_voltage_v:5.0f} V/turn".format(**row))
    for row in scan:
        if row["in_window"]:
            print("{gas} {fill_mtorr:>4.0f}mTorr {b_ext_t:.2f}T  T={total_t_ev:7.1f}eV s={s_bar:4.2f} "
                  "tau={tau_phi_classical_s:9.3e}s bias={bias_t:+.4f}T tq={quarter_period_s:.2e}s "
                  "E={bank_energy_j:7.0f}J V={bank_voltage_v:7.0f}V C={bank_capacitance_f:.2e}F "
                  "PI>={pi_te_min_ev:.2f}eV".format(**row))


if __name__ == "__main__":
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))
