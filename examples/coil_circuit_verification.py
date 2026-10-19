"""Coil circuit model against independent closed forms and published NUCTE-III hardware numbers."""
import json
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np
from plasma_bench.coil_circuit import Strap, section_coil, nagaoka_inductance, Module, Crowbar, simulate
from plasma_bench.frc_scaling import wheeler_inductance


def inductance_checks():
    """Uniform-current filament model -> Nagaoka's exact sheet value (second order in filament width);
    equal-flux (perfect conductor) inductance alongside."""
    rows = []
    for a, l in ((.1, .5), (.17, 1.4), (.08, .05), (.165, 1.2)):
        ref = nagaoka_inductance(a, l)
        for n in (16, 32, 64, 128):
            c = section_coil([Strap(a, 0., l)], n)
            rows.append(dict(radius_m=a, length_m=l, filaments=n, nagaoka_nh=ref*1e9,
                             uniform_relative_error=c.uniform_inductance()/ref-1,
                             equal_flux_over_nagaoka=c.inductance[0, 0]/ref))
    return rows


def multiturn_checks():
    """N series turns of width l/N: L -> N^2 x Nagaoka (uniform current per turn; equal flux only across
    each turn's own width, so the limit is approached from below as the turn width shrinks)."""
    out = []
    a, l = .1, .3
    for n_turns in (5, 10, 20, 40):
        w = l/n_turns
        straps = [Strap(a, -l/2+(k+.5)*w, w, turn=k) for k in range(n_turns)]
        coil = section_coil(straps, 8)
        out.append(dict(turns=n_turns, inductance_over_n2_nagaoka=float(coil.inductance[0, 0]/(n_turns**2*nagaoka_inductance(a, l)))))
    return out


def nucte_coil():
    """NUCTE-III coil (case-1 two-radius model, 28 straps x 5 cm) against the published 0.07 uH."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from examples.nucte_formation import Nucte, strap_layout
    out = {}
    for label, c in (("case1_two_radius", Nucte()), ("case6_graded", replace(Nucte(), coil_radius=.18, end_radius=.1267, taper_start=.1))):
        coil = section_coil([Strap(r, z, c.strap_width) for r, z in strap_layout(c)], 8)
        bz = sum(coil.field(0., 0.)[0][2] for _ in (0,))
        out[label] = dict(equal_flux_nh=coil.inductance[0, 0]*1e9, uniform_nh=coil.uniform_inductance()*1e9,
                          centre_field_t_per_ka=float(bz)*1e3)
    out["published_nh"] = 70.
    out["wheeler_r017_l14_nh"] = wheeler_inductance(.17, 1.4)*1e9
    return out


def lc_checks():
    """Single module on a single-section coil: lossless LC, underdamped RLC, ideal crowbar at V_C = 0."""
    coil = section_coil([Strap(.1, 0., .5)], 16)
    lc_ = coil.inductance[0, 0]
    c, v, ls, r = 2e-6, 10e3, 40e-9, 5e-3
    out = {}
    res = simulate(coil.inductance, [Module(c, v, ls)], 20e-6, samples=4001)
    l_tot = lc_+ls
    w = 1/np.sqrt(l_tot*c)
    exact = v*np.sqrt(c/l_tot)*np.sin(w*res.time)
    out["lc_max_relative_error"] = float(np.abs(res.section_current[:, 0]-exact).max()/np.abs(exact).max())
    out["lc_energy_residual_j"] = float(np.abs(res.energy_residual(coil.inductance)).max())
    res = simulate(coil.inductance, [Module(c, v, ls, resistance=r)], 20e-6, samples=4001)
    alpha = r/(2*l_tot); wd = np.sqrt(1/(l_tot*c)-alpha**2)
    exact = v/(wd*l_tot)*np.exp(-alpha*res.time)*np.sin(wd*res.time)
    out["rlc_max_relative_error"] = float(np.abs(res.section_current[:, 0]-exact).max()/np.abs(exact).max())
    out["rlc_energy_residual_relative"] = float(np.abs(res.energy_residual(coil.inductance)).max()/(.5*c*v*v))
    # Ideal crowbar (L_cb = R_cb = 0) at the capacitor-voltage zero, coil resistance R_c:
    # coil current decays as exp(-R_c t/L_c); the module rings in its own loop at 1/sqrt(L_s C).
    rc_ = 2e-4
    res = simulate(coil.inductance, [Module(c, v, ls, crowbar=Crowbar())], 60e-6, coil_resistance=rc_, samples=8001)
    tc = res.crowbar_times[0]
    after = res.time > tc
    tau = res.time[after]-tc
    i0 = np.interp(tc, res.time, res.section_current[:, 0])
    out["crowbar_time_s"] = tc
    out["crowbar_quarter_period_s"] = float(np.pi/2*np.sqrt(l_tot*c))
    out["crowbar_coil_decay_max_relative_error"] = float(np.abs(res.section_current[after, 0]/(i0*np.exp(-rc_*tau/lc_))-1).max())
    ws = 1/np.sqrt(ls*c)
    exact_m = i0*np.cos(ws*tau)
    out["crowbar_module_ring_max_relative_error"] = float(np.abs(res.module_current[after, 0]-exact_m).max()/i0)
    out["capacitor_reversal_v"] = float(res.capacitor_voltage[after, 0].min())
    out["capacitor_reversal_expected_v"] = float(-i0*np.sqrt(ls/c))
    out["crowbar_energy_residual_relative"] = float(np.abs(res.energy_residual(coil.inductance)).max()/(.5*c*v*v))
    # Unidirectional start switch: after the crowbar the module rings I_pk cos(w_s t) to its first zero,
    # a quarter period pi/2 sqrt(L_s C) later, then the capacitor is isolated holding -I_pk sqrt(L_s/C).
    res = simulate(coil.inductance, [Module(c, v, ls, crowbar=Crowbar(), unidirectional=True)], 30e-6,
                   coil_resistance=rc_, samples=6001)
    tc, to = res.crowbar_times[0], res.open_times[0]
    i0 = np.interp(tc, res.time, res.section_current[:, 0])
    out["switch_open_time_error_s"] = float(to-(tc+np.pi/2*np.sqrt(ls*c)))
    out["switch_held_voltage_relative_error"] = float(res.capacitor_voltage[-1, 0]/(-i0*np.sqrt(ls/c))-1)
    out["switch_energy_residual_relative"] = float(np.abs(res.energy_residual(coil.inductance)).max()/(.5*c*v*v))
    return out


def section_checks():
    """Mirror-image halves driven identically equal the common-collector coil; skewed and coupled
    sections still close energy."""
    straps = [Strap(.12, z, .05) for z in np.arange(-.575, .6, .05)]
    one = section_coil(straps, 8)
    two = section_coil([replace(s, section=int(s.z > 0)) for s in straps], 8)
    mod = dict(capacitance=1e-6, voltage=10e3, inductance=60e-9)
    res1 = simulate(one.inductance, [Module(**mod), Module(**mod)], 8e-6, samples=2001)
    res2 = simulate(two.inductance, [Module(**mod, section=0), Module(**mod, section=1)], 8e-6, samples=2001)
    total2 = res2.section_current.sum(axis=1)
    out = dict(common_collector_nh=one.inductance[0, 0]*1e9, halves_matrix_nh=(two.inductance*1e9).tolist(),
               halves_parallel_equivalent_nh=float(1e9/np.linalg.inv(two.inductance).sum()),
               symmetric_halves_vs_common_max_relative=float(np.abs(total2-res1.section_current[:, 0]).max()
                                                             / np.abs(res1.section_current).max()))
    skew = simulate(two.inductance, [Module(**mod, section=0, crowbar=Crowbar(resistance=1e-3)),
                                     Module(**mod, section=1, fire_time=100e-9, crowbar=Crowbar(resistance=1e-3))],
                    40e-6, coil_resistance=1e-4, samples=4001)
    out["skewed_energy_residual_relative"] = float(np.abs(skew.energy_residual(two.inductance)).max()/1e2)
    return out


def nucte_bank():
    """Published NUCTE-III circuit (Kuroyanagi 2014; Takahashi 2000) through the model.
    Strays are not published: the main and PI strays follow from the published periods, so those are
    consistency statements; the bias rise is an independent check."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from examples.nucte_formation import Nucte, strap_layout
    c = Nucte()
    coil = section_coil([Strap(r, z, c.strap_width) for r, z in strap_layout(c)], 8)
    lc_ = coil.inductance[0, 0]
    b_per_a = float(coil.field(0., 0.)[0][2])
    out = dict(coil_nh=lc_*1e9, centre_field_t_per_ka=b_per_a*1e3)
    l_main = (2*4e-6/np.pi)**2/67.5e-6
    out["main_total_nh_from_4us_rise"] = l_main*1e9
    out["main_stray_nh"] = (l_main-lc_)*1e9
    l_pi = 1/((2*np.pi*124e3)**2*11.25e-6)
    out["pi_total_nh_from_124khz"] = l_pi*1e9
    out["pi_stray_nh"] = (l_pi-lc_)*1e9
    for kv in (32., 20., 13.):
        res = simulate(coil.inductance, [Module(67.5e-6, kv*1e3, l_main-lc_, crowbar=Crowbar())], 6e-6, samples=3001)
        out[f"main_{kv:g}kV_peak_centre_field_t"] = float(res.section_current[:, 0].max()*b_per_a)
    bias = simulate(coil.inductance, [Module(2000e-6, 1., 2.4e-6, resistance=8e-3)], 300e-6, samples=6001)
    k = int(np.argmax(bias.section_current[:, 0]))
    out["bias_time_to_peak_us"] = float(bias.time[k]*1e6)
    out["bias_published_rise_us"] = 90.
    out["bias_kv_per_0p09T"] = float(.09/(bias.section_current[k, 0]*b_per_a)/1e3)
    out["bias_peak_current_ka_at_0p09T"] = float(.09/b_per_a/1e3)
    return out


def main(out=Path("outputs/coil_circuit_verification")):
    out.mkdir(parents=True, exist_ok=True)
    result = dict(inductance=inductance_checks(), multiturn=multiturn_checks(), nucte_coil=nucte_coil(), lc=lc_checks(),
                  sections=section_checks(), nucte_bank=nucte_bank())
    (out/"results.json").write_text(json.dumps(result, indent=1))
    for r in result["inductance"]:
        print("a=%.3f l=%.2f n=%3d Nagaoka %.2f nH  uniform err %+.2e  equal-flux/Nagaoka %.4f" % (
            r["radius_m"], r["length_m"], r["filaments"], r["nagaoka_nh"], r["uniform_relative_error"], r["equal_flux_over_nagaoka"]))
    for key in ("multiturn", "nucte_coil", "lc", "sections", "nucte_bank"):
        print(key, json.dumps(result[key], indent=1))


if __name__ == "__main__":
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))
