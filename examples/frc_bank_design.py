"""Bank, crowbar, timing and mirror numbers for the 12-inch design through the coil circuit model.

Everything here is vacuum circuit + coil-field arithmetic (plasma loading is estimated separately);
the formation consequences come from examples/nucte_formation.py runs."""
import json
import sys
from pathlib import Path
import numpy as np
from plasma_bench.coil_circuit import Strap, section_coil, Module, Crowbar, simulate

R_COIL, LENGTH, WIDTH = .165, 1.2, .05        # single-turn coil hugging a 12" ID tube (FRC_DESIGN_WINDOW.md)
MODULES, C_MOD, L_MOD = 8, 6e-6, 30e-9        # baseline bank: eight ~6 uF modules, 30 nH stray each
BIAS = dict(capacitance=2000e-6, inductance=2.4e-6, resistance=8e-3)   # NUCTE-III-like slow bias bank


def straps(radius=lambda z: R_COIL, sections=1, mirror=None):
    """Single-turn main coil (24 x 5 cm straps) on `sections` collectors; optional mirror pair (one series
    circuit, its own section): dict(turns, radius, width, gap) placed just beyond both coil ends."""
    zs = -LENGTH/2+WIDTH*(np.arange(int(round(LENGTH/WIDTH)))+.5)
    per = int(np.ceil(zs.size/sections))
    out = [Strap(float(radius(z)), float(z), WIDTH, section=k//per) for k, z in enumerate(zs)]
    if mirror:
        sec = sections
        turn = 0
        for side in (-1, 1):
            for j in range(mirror["turns"]):
                z = side*(LENGTH/2+mirror["gap"]+(j+.5)*mirror["width"])
                out.append(Strap(mirror["radius"], float(z), mirror["width"], section=sec, turn=turn))
                turn += 1
    return out


def profile(coil, currents, z=np.linspace(0, 1.0, 1001)):
    """On-axis Bz(z) for given section currents, mirror ratio and NUCTE field index (z/B dB/dz averaged
    between the mirror points)."""
    fields = coil.field(np.zeros_like(z), z)
    b = sum(i*f[2] for i, f in zip(currents, fields))
    k = int(np.argmax(b))
    nz = z[:k+1]*np.gradient(b, z)[:k+1]/b[:k+1]
    return dict(b0=float(b[0]), mirror_ratio=float(b[k]/b[0]), mirror_point_m=float(z[k]),
                field_index=float(np.trapezoid(nz, z[:k+1])/z[k]) if k > 0 else 0.)


def main_bank(coil, b_peak, b_bias, sections=None, skew=None, l_cb=20e-9, hold=300e-6, unidirectional=True):
    """Bias + main banks for a target vacuum centre field (bias at t=0, peak after the rise).
    Charge voltages from linear superposition (one correction pass). Crowbar resistance per module for an
    L/R hold `hold` of the coil + crowbar loop."""
    n_s = coil.sections
    bz = np.array([f[2] for f in coil.field(0., 0.)], dtype=float)
    main_sec = list(range(n_s)) if sections is None else sections
    per = MODULES//len(main_sec)
    skew = skew or {}
    l_coil = 1/np.linalg.inv(coil.inductance[np.ix_(main_sec, main_sec)]).sum()
    l_cb_tot = l_cb/MODULES
    r_cb = (l_coil+l_cb_tot)/hold*MODULES

    def banks(v_main, v_bias, lead):
        mods = [Module(BIAS["capacitance"]/len(main_sec), v_bias, BIAS["inductance"]*len(main_sec),
                       BIAS["resistance"]*len(main_sec), section=s, fire_time=-lead, name=f"bias{s}") for s in main_sec]
        k = 0
        for s in main_sec:
            for _ in range(per):
                mods.append(Module(C_MOD, v_main, L_MOD, section=s, fire_time=skew.get(k, 0.),
                                   crowbar=Crowbar(l_cb, r_cb), unidirectional=unidirectional, name=f"main{k}"))
                k += 1
        return mods

    centre = lambda res: res.section_current@bz
    probe = simulate(coil.inductance, [m for m in banks(0., 1., 0.) if m.name.startswith("bias")], 400e-6, samples=4001)
    lead = float(probe.time[np.argmax(np.abs(probe.section_current).sum(axis=1))])
    v_bias = b_bias/float(np.interp(lead, probe.time, centre(probe))) if b_bias else 0.
    unit = simulate(coil.inductance, [m for m in banks(1., 0., 0.) if m.name.startswith("main")], 8e-6, samples=4001)
    v_main = (b_peak-b_bias)/float(centre(unit).max())
    grid = lambda end: np.concatenate((np.linspace(-lead, 0, 201)[:-1], np.linspace(0, end, 4001)))
    for _ in range(2):
        res = simulate(coil.inductance, banks(v_main, v_bias, lead), 30e-6, samples=grid(30e-6))
        f0 = float(np.interp(0., res.time, centre(res)))
        v_main *= (b_peak-f0)/(float(centre(res)[res.time >= 0].max())-f0)
    res = simulate(coil.inductance, banks(v_main, v_bias, lead), 30e-6, samples=grid(30e-6))
    return res, dict(v_main=v_main, v_bias=v_bias, lead=lead, r_cb=r_cb, l_coil=l_coil, bz=bz)


def bank_summary(res, info):
    mains = [k for k, m in enumerate(res.modules) if m.name.startswith("main")]
    t = res.time
    cf = res.section_current@info["bz"]
    tp = res.crowbar_times[mains[0]]
    after = t >= tp
    i_mod = res.module_current[:, mains]
    di = np.gradient(i_mod, t, axis=0)
    hold = lambda dt: float(np.interp(tp+dt, t, cf)/np.interp(tp, t, cf))
    return dict(charge_kv=info["v_main"]/1e3, bias_kv=info["v_bias"]/1e3, bank_energy_kj=sum(.5*res.modules[k].capacitance*res.modules[k].voltage**2 for k in mains)/1e3,
                coil_nh=info["l_coil"]*1e9, crowbar_s=tp, peak_module_ka=float(np.abs(i_mod).max()/1e3),
                peak_module_di_dt_ka_per_us=float(np.abs(di[t < tp]).max()/1e9),
                capacitor_reversal_kv=float(res.capacitor_voltage[:, mains].min()/1e3),
                reversal_fraction=float(-res.capacitor_voltage[:, mains].min()/info["v_main"]),
                field_at_crowbar_t=float(np.interp(tp, t, cf)), held_fraction_2us=hold(2e-6), held_fraction_10us=hold(10e-6),
                held_fraction_25us=hold(min(25e-6, t[-1]-tp)), crowbar_r_mohm_per_module=info["r_cb"]*1e3)


def skew_common(coil, b_peak, b_bias, late):
    """Common collector: one module fires `late` after the other seven; total field and switch currents."""
    ref, info = main_bank(coil, b_peak, b_bias)
    res, _ = main_bank(coil, b_peak, b_bias, skew={0: late})
    cf0, cf = ref.section_current@info["bz"], res.section_current@info["bz"]
    t = ref.time
    early = (t > 0) & (t < 1.5e-6)
    mains = [k for k, m in enumerate(res.modules) if m.name.startswith("main")]
    return dict(late_ns=late*1e9, max_field_change_t=float(np.abs(np.interp(t, res.time, cf)-cf0)[early].max()),
                relative_to_field_at_0p5us=float(np.abs(np.interp(t, res.time, cf)-cf0)[early].max()/abs(np.interp(.5e-6, t, cf0))),
                peak_current_late_ka=float(np.abs(res.module_current[:, mains[0]]).max()/1e3),
                peak_current_others_ka=float(np.abs(res.module_current[:, mains[1]]).max()/1e3))


def skew_sections(b_peak, b_bias, late, n_sections=8):
    """Separate sections (one module each): modules on z > 0 fire `late` later. Axial asymmetry of the
    vacuum field on the wall radius during the rise."""
    coil = section_coil(straps(sections=n_sections), 8)
    res, info = main_bank(coil, b_peak, b_bias, skew={k: late for k in range(n_sections//2, n_sections)})
    z = np.linspace(-.55, .55, 111)
    fields = coil.field(np.full_like(z, .14), z)
    out = {}
    for t in (.2e-6, .5e-6, 1e-6, 2e-6):
        cur = [np.interp(t, res.time, res.section_current[:, s]) for s in range(coil.sections)]
        b = sum(i*f[2] for i, f in zip(cur, fields))
        out[f"{t*1e6:g}us"] = float(np.abs(b-b[::-1]).max()/np.abs(b).max())
    return dict(late_ns=late*1e9, left_right_asymmetry=out)


def main(out=Path("outputs/frc_bank_design")):
    out.mkdir(parents=True, exist_ok=True)
    result = {}
    base = section_coil(straps(), 8)
    result["coil"] = dict(inductance_nh=float(base.inductance[0, 0]*1e9), centre_field_mt_per_ka=float(base.field(0., 0.)[0][2]*1e6),
                          **profile(base, [1.]))
    result["banks"] = {}
    for b_peak, b_bias in ((.2, -.05), (.3, -.05), (.3, -.07)):
        for l_cb in (8e-9, 20e-9, 40e-9):
            res, info = main_bank(base, b_peak, b_bias, l_cb=l_cb)
            result["banks"][f"B{b_peak}_bias{b_bias}_Lcb{l_cb*1e9:g}nH"] = bank_summary(res, info)
    result["skew_common"] = [skew_common(base, .3, -.05, late) for late in (50e-9, 200e-9, 500e-9)]
    result["skew_sections"] = [skew_sections(.3, -.05, late) for late in (10e-9, 50e-9, 200e-9)]
    mirrors = {}
    for r0 in (.20, .22, .24):
        coil = section_coil(straps(radius=lambda z, r0=r0: r0+(R_COIL-r0)*np.clip((abs(z)-.1)/(LENGTH/2-.1), 0, 1)), 8)
        res, info = main_bank(coil, .3, -.05)
        mirrors[f"graded_{r0:.2f}_to_{R_COIL}"] = dict(**profile(coil, [1.]), inductance_nh=float(coil.inductance[0, 0]*1e9),
                                                       **{k: v for k, v in bank_summary(res, info).items() if k in ("charge_kv", "bank_energy_kj", "peak_module_ka")})
    for turns in (2, 4):
        coil = section_coil(straps(mirror=dict(turns=turns, radius=.17, width=.02, gap=.01)), 8)
        # Mirror pair on its own module, fired with the main bank: scan its share for R_M.
        bz = np.array([f[2] for f in coil.field(0., 0.)])
        for ratio in (.2, .4, .6):
            prof = profile(coil, [1., ratio/turns])
            mirrors[f"mirror_pair_{turns}turn_I{ratio:.1f}"] = dict(**prof, mirror_current_fraction_of_main=ratio/turns,
                                                                    mirror_inductance_nh=float(coil.inductance[1, 1]*1e9),
                                                                    mutual_nh=float(coil.inductance[0, 1]*1e9))
    result["mirrors"] = mirrors
    (out/"results.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))
    return result


if __name__ == "__main__":
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))
