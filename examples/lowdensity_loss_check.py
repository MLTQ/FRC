"""Low-density heat sinks of the 3D hybrid (frc-ibl.11.9): what the numerical removal stands in for.

`python -m examples.lowdensity_loss_check RUN [t_us ...]` reads saved snapshots of a completed hybrid3d_frc run
(snap_NNNN.npz: b_response, n, pe) and its setup, and writes RUN/lowdensity_loss_check.json. Post-processing
only: no solver state changes and no new runs.

For every gas node it rebuilds the solver's own Joule heating, q = edge average of eta J^2. Here J is
curl(B_response)/mu0 on the plasma edges (touching gas, not on a conductor), and eta is the solver's NRL Spitzer
1.03e-4 lnL max(Te, 1 eV)^-1.5 with Te = pe/(n e). It then splits the nodes into plasma (n >= floor),
thin (0 < n < floor) and empty (n <= 0, with gas). Three questions:

1. **Where is heat removed?** The removed part of q is q (1 - (n/max(n, floor))^gamma), as in EnergyElectrons.
2. **Are there carriers for that current?** The macro-ions per node are n f V/w (w = mean macro-ion weight).
   The electron drift is v_d = |J|/(n e) against the thermal speed v_te = sqrt(e Te/m_e). Where v_d > v_te, or
   there is less than one macro-ion, the current has no resolved carriers. Its Joule heat is then a property of
   the floor, not of a plasma.
3. **What would physics dissipate instead?** If the heat stayed in the electrons, Te would rise, Spitzer eta
   would fall as Te^-1.5, and the heat would leave along open field lines. Holding J fixed, solve per node
       eta(T) J^2 = (3/2) n T f v_te(T)/L,   f = 0.1 (flux-limited parallel conduction)
   with L the axial distance to the nearer tube end. Open edge lines run axially to the ends; this is declared
   and not traced. The steady dissipation eta(T_ss) J^2 is the heat a self-heating, conducting thin plasma
   would actually remove. It is compared with the numerical removal. Frozen J makes this a bound on rates,
   not a replacement model; the counterfactual changes no state.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, m_e as ME
from plasma_bench.hybrid3d import curl_faces, edges_to_nodes
from plasma_bench.hybrid3d_device import Device, device_box, cell_types, masks, GAS

LNL, GAMMA, F_COND = 10., 5/3, .1


def eta_spitzer(te):
    return 1.03e-4*LNL*np.maximum(te, 1.)**-1.5


def steady_temperature(j2, n, length):
    """T [eV] with eta(T) j2 = 1.5 n e T f v_te(T)/L (bisection in log T, 1 eV .. 100 keV), per node."""
    lo, hi = np.zeros_like(n), np.full_like(n, np.log(1e5))

    def balance(lnT):
        t = np.exp(lnT)
        return eta_spitzer(t)*j2-1.5*n*QE*t*F_COND*np.sqrt(QE*t/ME)/length
    for _ in range(60):
        mid = .5*(lo+hi)
        pos = balance(mid) > 0         # heating exceeds loss: hotter
        lo, hi = np.where(pos, mid, lo), np.where(pos, hi, mid)
    return np.exp(.5*(lo+hi))


def analyse(run, t_us):
    run = Path(run)
    info = json.loads((run/"setup.json").read_text())
    dv = info["device"]
    dev = Device(dv["wall"], dv["half"], dv["coil"], dv["end_gap"], dv.get("inner", 0.))
    box = device_box(dev, info["dx"], info["dz"])
    if list(box.shape) != list(info["box"]):
        raise ValueError(f"rebuilt box {box.shape} != saved {info['box']}")
    types = cell_types(box, dev)
    _, _, edge_mask = masks(types)
    gas = (types == GAS).astype(float)
    frac = gas.copy()
    for a in range(3):
        frac = .5*(frac+np.roll(frac, 1, axis=a))
    live = frac > 0
    vol = box.cell_volume
    floor = info["floor_m3"]
    hist = json.loads((run/"history.json").read_text())
    w_mean = hist[0]["particles"]/info["particles"]
    z = box.lo[2]+box.d[2]*np.arange(box.shape[2])
    length = np.maximum(dev.tube_half_length-np.abs(z), box.d[2])[None, None, :]*np.ones(box.shape)
    step = info["snapshot_interval_s"]
    rows = []
    for t in t_us:
        path = run/f"snap_{int(round(t*1e-6/step)):04d}.npz"
        s = np.load(path)
        br = [s["b_response"][a].astype(float) for a in range(3)]
        n = s["n"].astype(float); pe = s["pe"].astype(float)
        te = pe/(np.maximum(n, 1.)*QE)
        j_edges = [c/MU0*edge_mask[a] for a, c in enumerate(curl_faces(*br, box.d))]
        eta = eta_spitzer(te)*live
        eta_e = [.5*(eta+np.roll(eta, -1, axis=a))*edge_mask[a] for a in range(3)]
        q = np.sum(edges_to_nodes(*[eta_e[a]*j_edges[a]**2 for a in range(3)]), axis=0)
        jn = edges_to_nodes(*j_edges)
        jn = np.where(live, jn/np.where(live, frac, 1.), 0.)
        j2 = np.sum(jn**2, axis=0)
        retained = q*(np.maximum(n, 0.)/np.maximum(np.maximum(n, floor), 1.))**GAMMA
        removed = q-retained
        cls = {"plasma": live & (n >= floor), "thin": live & (n > 0) & (n < floor), "empty": live & (n <= 0)}
        carriers = n*frac*vol/w_mean
        vte = np.sqrt(QE*np.maximum(te, 1.)/ME)
        vd = np.sqrt(j2)/(np.maximum(n, 1.)*QE)
        unresolved = live & ((carriers < 1) | (vd > vte))
        ts = steady_temperature(j2, np.maximum(n, 1.), length)
        q_phys = np.where(cls["thin"] | cls["plasma"], eta_spitzer(ts)*j2, 0.)
        row = dict(time_us=float(s["time_s"])*1e6)
        for name, m in cls.items():
            row[name] = dict(nodes=int(m.sum()), joule_w=float(q[m].sum()*vol), removed_w=float(removed[m].sum()*vol),
                             mean_macro_ions=float(carriers[m].mean()) if m.any() else 0.,
                             removed_from_unresolved_w=float(removed[m & unresolved].sum()*vol))
        thin = cls["thin"]
        row["thin_counterfactual"] = dict(
            removed_w=float(removed[thin].sum()*vol),
            self_heating_dissipation_w=float(q_phys[thin].sum()*vol),
            median_te_now_ev=float(np.median(te[thin])) if thin.any() else None,
            median_te_steady_ev=float(np.median(ts[thin])) if thin.any() else None,
            median_vd_over_vte=float(np.median((vd/vte)[thin])) if thin.any() else None)
        rows.append(row)
        print(json.dumps(row), flush=True)
    out = dict(run=str(run), floor_m3=floor, mean_macro_ion_weight=w_mean, flux_limiter=F_COND, coulomb_log=LNL, rows=rows)
    (run/"lowdensity_loss_check.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    analyse(sys.argv[1], [float(v) for v in sys.argv[2:]] or [.2, .5, 1., 2., 4., 10.])
