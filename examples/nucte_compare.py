"""Post-process a nucte_formation run into the quantities NUCTE-III published (Takahashi 2000)."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, atomic_mass as AMU

D_MASS = 2.01410*AMU
# Takahashi et al. 2000, Table 1 (experiments on six coil configurations) and text; units SI.
# tau_ax: onset of axial contraction; tau_Max: contraction peak at the midplane; tau_eq: equilibrium.
# n is the interferometer line average over the separatrix, int(n dl)/(2 r_s(0)); T = Te+Ti is from
# pressure balance with <beta> = 1 - x_s^2/2 (x_s = r_s/r_coil). Lifetime is limited by the n=2 mode.
NUCTE_CASE1 = dict(separatrix_mid_m=.062, separatrix_length_m=.62, contraction_onset_s=6.5e-6, axial_speed_m_s=3.8e4,
                   contraction_peak_s=12.5e-6, equilibrium_s=30e-6, lifetime_s=65e-6, n2_onset_s=28e-6,
                   density_m3=2.2e21, total_temperature_ev=120., mirror_ratio=1.10, field_index=.03)
# Case 6 text: reconnection complete by 3 us, contraction front at the midplane at 6 us, equilibrium
# at 10 us with r_s ~ 8 cm and l_s ~ 40 cm.
NUCTE_CASE6 = dict(separatrix_mid_m=.080, separatrix_length_m=.38, contraction_onset_s=3.0e-6, axial_speed_m_s=8.5e4,
                   contraction_peak_s=7.0e-6, equilibrium_s=12e-6, lifetime_s=50e-6, n2_onset_s=22e-6,
                   density_m3=2.2e21, total_temperature_ev=120., mirror_ratio=1.80, field_index=.30)


def reference(setup):
    """Case 6 for strong-mirror coils (end straps well inside the centre straps), else case 1."""
    case = setup["case"]
    return ("case 6", NUCTE_CASE6) if (case["coil_radius"]/case["end_radius"])**2 > 1.4 else ("case 1", NUCTE_CASE1)


def coil_profile(case, points=1001):
    """Axial vacuum field of the modeled coil: mirror ratio B_max/B(0) and NUCTE's field index
    n_z = (z/B) dB/dz averaged between the mirror points (Takahashi 2000)."""
    from dataclasses import fields
    from examples.nucte_formation import Nucte, equal_flux_currents
    from plasma_bench.electromagnetics import loop_field
    c = Nucte(**{f.name: case[f.name] for f in fields(Nucte) if f.name in case})
    fil, cur = equal_flux_currents(c)
    z = np.linspace(0, c.half_length, points)
    b = np.zeros_like(z)
    for f, i in zip(fil, cur):
        for rr, zz in f:
            b += i*loop_field(np.zeros_like(z), z, rr, zz)[2]/len(f)
    k = int(np.argmax(b))
    nz = z[:k+1]*np.gradient(b, z)[:k+1]/b[:k+1]
    return dict(mirror_ratio=float(b[k]/b[0]), mirror_point_m=float(z[k]),
                field_index=float(np.trapezoid(nz, z[:k+1])/z[k]))


def load(run):
    run = Path(run)
    setup = json.loads((run/"setup.json").read_text())
    data = np.load(run/"snapshots.npz")
    rf, zf = data["rf"], data["zf"]
    nr, nz = len(rf)-1, len(zf)-1
    frames = []
    for key in sorted(k for k in data.files if k.startswith("t")):
        v = data[key]
        n_psi = (nr+1)*(nz+1)
        psi = v[:n_psi].reshape(nr+1, nz+1)
        rho = v[n_psi:n_psi+nr*nz].reshape(nr, nz)
        p = v[n_psi+nr*nz:].reshape(nr, nz)
        frames.append((int(key[1:]), psi, rho, p))
    history = json.loads((run/"history.json").read_text())
    return setup, rf, zf, frames, history


def separatrix_profile(rf, psi):
    """r_s at every axial node: zero of total flux outward from a reversed axis; nan if unreversed."""
    rs = np.full(psi.shape[1], np.nan)
    for j in range(psi.shape[1]):
        col = psi[:, j]
        if col[1] < 0:
            k = np.nonzero(col[1:] >= 0)[0]
            if k.size:
                i = k[0]+1
                rs[j] = rf[i-1]+(rf[i]-rf[i-1])*(-col[i-1])/(col[i]-col[i-1])
            else:
                rs[j] = np.inf
    return rs


def closed_region(rf, zf, psi):
    """Closed poloidal flux: reversed-axis segment whose separatrix stays off the wall, bounded by unreversed ends."""
    rs = separatrix_profile(rf, psi)
    reversed_ = np.isfinite(rs) & (rs < rf[-1])
    if not reversed_.any():
        return None
    mid = len(zf)//2
    # contiguous reversed segment containing (or nearest to) the midplane
    idx = np.nonzero(reversed_)[0]
    centre = idx[np.argmin(abs(idx-mid))]
    lo = hi = centre
    while lo > 0 and reversed_[lo-1]:
        lo -= 1
    while hi < len(zf)-1 and reversed_[hi+1]:
        hi += 1
    ends_open = (lo == 0) or (hi == len(zf)-1) or np.isinf(rs[max(lo-1, 0)]) or np.isinf(rs[min(hi+1, len(zf)-1)])
    o_point = float(-psi[:, lo:hi+1].min())
    return dict(lo=int(lo), hi=int(hi), length_m=float(zf[hi]-zf[lo]), closed=not ends_open,
                separatrix_mid_m=float(rs[mid]) if np.isfinite(rs[mid]) else None,
                max_separatrix_m=float(np.nanmax(np.where(np.isfinite(rs[lo:hi+1]), rs[lo:hi+1], np.nan))),
                trapped_flux_wb=2*np.pi*o_point, profile=rs)


def plasma_inside(rf, zf, rho, p, rs, lo, hi):
    """Volume-average density and total temperature (Te+Ti = p/(n e)) inside r < r_s(z)."""
    rc = .5*(rf[:-1]+rf[1:])
    dz = zf[1]-zf[0]
    vol = np.pi*np.diff(rf**2)[:, None]*np.full(len(zf)-1, dz)[None, :]
    inside = np.zeros_like(rho, dtype=bool)
    for j in range(lo, hi):
        r_s = np.nanmin([rs[j], rs[j+1]])
        inside[:, j] = rc < r_s
    if not inside.any():
        return None, None
    n = rho/D_MASS
    nbar = float(np.sum(n*vol*inside)/np.sum(vol*inside))
    temperature = float(np.sum(p*vol*inside)/np.sum(n*vol*inside)/QE)  # p = n_i k(Ti+Te), Z=1
    return nbar, temperature


def chord_density(rf, rho, j):
    """Side-on chord through the axis at axial cell j: integral n dl = 2 * integral_0^rw n dr."""
    return float(2*np.sum(rho[:, j]/D_MASS*np.diff(rf)))


def summarize(run):
    setup, rf, zf, frames, history = load(run)
    rows = []
    for k, psi, rho, p in frames:
        t = history[k]["time_s"] if k < len(history) else np.nan
        region = closed_region(rf, zf, psi)
        mid_cell = rho.shape[1]//2
        row = dict(time_s=t, chord_density_m2=chord_density(rf, rho, mid_cell))
        if region:
            nbar, temp = plasma_inside(rf, zf, rho, p, region["profile"], region["lo"], region["hi"])
            row.update({k2: v for k2, v in region.items() if k2 != "profile"}, density_m3=nbar, total_temperature_ev=temp)
            row["profile"] = [None if not np.isfinite(x) else float(x) for x in region["profile"]]
        rows.append(row)
    return setup, rf, zf, rows


def decay(rows, start_s=8e-6):
    """Trapped-flux decay after formation: least-squares e-folding time of ln(phi) over closed frames
    with t >= start_s, the last time a closed FRC exists, and the flux there. None if too few frames."""
    good = [(r["time_s"], r["trapped_flux_wb"]) for r in rows
            if r.get("closed") and r.get("trapped_flux_wb", 0) > 0 and r["time_s"] >= start_s]
    last = max((r["time_s"] for r in rows if r.get("closed")), default=None)
    if len(good) < 3:
        return dict(flux_efold_s=None, last_closed_s=last)
    t, phi = np.array(good).T
    slope = np.polyfit(t, np.log(phi), 1)[0]
    return dict(flux_efold_s=float(-1/slope) if slope < 0 else None, last_closed_s=last,
                fit_window_s=[float(t[0]), float(t[-1])])


def main(run, out=None):
    setup, rf, zf, rows = summarize(run)
    for r in rows:
        print("t=%5.1f us closed=%s rs(0)=%s length=%s phi=%s n=%s T=%s chord=%.3g" % (
            r["time_s"]*1e6, r.get("closed"), None if r.get("separatrix_mid_m") is None else round(r["separatrix_mid_m"]*100, 2),
            None if "length_m" not in r else round(r["length_m"]*100, 1),
            None if "trapped_flux_wb" not in r else round(r["trapped_flux_wb"]*1e3, 3),
            None if r.get("density_m3") is None else "%.2e" % r["density_m3"],
            None if r.get("total_temperature_ev") is None else round(r["total_temperature_ev"], 1), r["chord_density_m2"]))
    label, ref = reference(setup)
    lifetime = decay(rows)
    coil = coil_profile(setup["case"])
    print("flux decay:", lifetime)
    print("coil: mirror ratio %.2f (published %.2f), field index %.3f (published %.2f)" % (
        coil["mirror_ratio"], ref["mirror_ratio"], coil["field_index"], ref["field_index"]))
    (Path(out or run)/"comparison.json").write_text(json.dumps(dict(reference=label, targets=ref, decay=lifetime, coil=coil,
                                                                    rows=rows), indent=1))
    return rows


if __name__ == "__main__":
    main(sys.argv[1], *(sys.argv[2:3]))
