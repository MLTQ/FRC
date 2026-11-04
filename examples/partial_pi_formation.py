"""Formation from a partially pre-ionized fill with a static neutral reservoir (frc-ibl.4 step 3; see .md).

    python -m examples.partial_pi_formation check
    python -m examples.partial_pi_formation run OUT alpha=0.1 cost_ev=30 [cx=3e-19 end_us=12 preset=frc12_uniform_h2]

`check` verifies apply_ionization against the exact logistic solution of dn_e/dt = K n_e n_n (n_e + n_n fixed,
constant K): mass and energy bookkeeping exact, first-order convergence in dt, and momentum and total energy
conserved apart from the declared cost for a moving cell. It also checks charge exchange alone (constant K_cx):
neutrals n_n0 exp(-n_i K t) and momentum p0 exp(-(n_n0/n_i)(1 - exp(-n_i K t))), with the total energy change
equal to the booked exchange loss.

`run` is examples/nucte_formation.main with the preset unchanged except for the fill split: the plasma starts as the
fraction alpha of the fill under the coil (fill_pa scaled by alpha, floor_fraction by 1/alpha so the declared floor
and the residual plasma beyond the coil keep their absolute densities), and the rest of the fill under the coil is a
static neutral reservoir, (1 - alpha) n_fill times the same PI profile. alpha = 1 leaves the reservoir empty and must
reproduce the preset run exactly.
"""
import json
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.constants import e as QE
import examples.nucte_formation as nf
from plasma_bench.neutral_reservoir import NeutralReservoirMHD, apply_ionization


def check():
    gamma, mass, k, cost = 5/3, 1.67e-27, 2e-14, 30.
    n_total, n0 = 2.2e20, 2.2e18
    mesh = SimpleNamespace(volumes=np.ones((1, 1)), rf=np.array([0., 1.]))
    te0 = 20.
    end = 3/(k*n_total)
    exact = n_total/(1+(n_total/n0-1)*np.exp(-k*n_total*end))
    rows = []
    for steps in (100, 200, 400, 800):
        dt = end/steps
        u = np.zeros((1, 1, 4))
        u[..., 0] = n0*mass
        u[..., 3] = 2*n0*te0*QE/(gamma-1)*50          # ample thermal energy: the cost cap never binds
        neutral = np.full((1, 1), n_total-n0)
        e0 = u[..., 3].sum()
        for _ in range(steps):
            u, _, neutral, booked = apply_ionization(mesh, u, None, np.zeros((1, 1)), np.zeros((1, 1), bool), neutral,
                                                     dt, mass, gamma, lambda t: np.full_like(t, k), cost)
        ne = float(u[..., 0].sum()/mass)
        rows.append(dict(steps=steps, n_e=ne, relative_error=abs(ne-exact)/exact,
                         mass_conservation=float((u[..., 0].sum()/mass+neutral.sum()-n_total)/n_total),
                         energy_minus_cost=float((e0-u[..., 3].sum())/((ne-n0)*cost*QE)-1)))
    order = [float(np.log2(a['relative_error']/b['relative_error'])) for a, b in zip(rows, rows[1:])]
    # moving cell: momentum fixed, kinetic drop to heat, total energy changes only by the cost
    u = np.zeros((1, 1, 4)); u[..., 0] = n0*mass; u[..., 1] = n0*mass*1e5; u[..., 3] = .5*n0*mass*1e10+1e3
    moved, _, _, booked = apply_ionization(mesh, u, None, np.zeros((1, 1)), np.zeros((1, 1), bool),
                                           np.full((1, 1), n0), 1e-6, mass, gamma, lambda t: np.full_like(t, k), cost)
    moving = dict(momentum_change=float(moved[..., 1].sum()-u[..., 1].sum()),
                  energy_change_plus_cost=float(moved[..., 3].sum()-u[..., 3].sum()+booked['cost_j']),
                  loading_heat_j=booked['loading_heat_j'])
    # charge exchange alone: n_i fixed, neutrals and momentum decay analytically
    n_i, n_n0, kcx, end_cx = 2e20, 1e20, 3e-14, 2/(2e20*3e-14)
    exact_neutral = n_n0*np.exp(-n_i*kcx*end_cx)
    exact_momentum = np.exp(-(n_n0/n_i)*(1-np.exp(-n_i*kcx*end_cx)))
    cx_rows = []
    for steps in (100, 200, 400, 800):
        dt = end_cx/steps
        u = np.zeros((1, 1, 4)); u[..., 0] = n_i*mass; u[..., 1] = n_i*mass*1e5; u[..., 3] = .5*n_i*mass*1e10+3*n_i*QE*50
        neutral = np.full((1, 1), n_n0); e0, p0, lost = u[..., 3].sum(), u[..., 1].sum(), 0.
        for _ in range(steps):
            u, _, neutral, booked = apply_ionization(mesh, u, None, np.zeros((1, 1)), np.zeros((1, 1), bool), neutral, dt,
                                                     mass, gamma, lambda t: np.zeros_like(t), 0.,
                                                     lambda t, v2: np.full_like(t, kcx))
            lost += booked['cx_energy_j']
        cx_rows.append(dict(steps=steps, neutral_error=float(abs(neutral.sum()/exact_neutral-1)),
                            momentum_error=float(abs(u[..., 1].sum()/p0/exact_momentum-1)),
                            mass_change=float(u[..., 0].sum()/(n_i*mass)-1),
                            energy_change_plus_booked=float((u[..., 3].sum()-e0+lost)/e0)))
    cx_order = [float(np.log2(a['momentum_error']/b['momentum_error'])) for a, b in zip(cx_rows, cx_rows[1:])]
    out = dict(exact_n_e=float(exact), convergence=rows, observed_order=order, moving_cell=moving,
               charge_exchange=cx_rows, charge_exchange_momentum_order=cx_order)
    print(json.dumps(out, indent=1))
    return out


def run(out, alpha=.1, cost_ev=30., end_us=12., preset='frc12_uniform_h2', cx=0.):
    from examples.experiments import PRESETS
    alpha = float(alpha)
    overrides = dict(PRESETS[preset])
    case = nf.replace(nf.Nucte(), **overrides)
    n_fill = 2*case.fill_pa/(nf.KB*case.gas_k)
    overrides.update(fill_pa=case.fill_pa*alpha, floor_fraction=case.floor_fraction/alpha, end_s=float(end_us)*1e-6)
    build, observe = nf.build, nf.observe

    def build_with_reservoir(c):
        model, setup = build(c)
        model.__class__ = NeutralReservoirMHD
        zc = .5*(model.mesh.zf[:-1]+model.mesh.zf[1:])
        shape = .5*(1-np.tanh((np.abs(zc)-c.pi_half_length)/c.pi_edge_width)) if c.pi_half_length > 0 else np.ones_like(zc)
        neutral = (1-alpha)*n_fill*np.ones((c.nr, 1))*shape[None, :]
        model.attach_reservoir(neutral, model.ion_mass, float(cost_ev), cx_cross_section=float(cx))
        setup.update(partial_pi=dict(alpha=alpha, cost_ev=float(cost_ev), cx_cross_section_m2=float(cx),
                                     cx_model='escape limit: each exchanged ion leaves as a fast neutral' if float(cx) > 0 else None,
                                     true_fill_ions_m3=n_fill,
                                     neutral_atoms=model.neutral_initial,
                                     rate='NIST-BEB H2 total ionization (Maxwellian), applied per atom-equivalent'))
        return model, setup

    def observe_with_reservoir(model):
        row = observe(model)
        if hasattr(model, 'reservoir_diagnostics'):
            row['reservoir'] = model.reservoir_diagnostics()
        return row
    nf.build, nf.observe = build_with_reservoir, observe_with_reservoir
    return nf.main(Path(out), **overrides)


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        check()
    else:
        args = dict(a.split('=', 1) for a in sys.argv[3:])
        run(sys.argv[2], **args)
