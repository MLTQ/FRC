"""0D pre-ionization energy and ionization budget for the 12-inch design (frc-ibl.4, step 1).

    python -m examples.pi_budget [OUT]        (default outputs/pi_budget; prints the tables)

What is exact and what is declared:
- Inventory: atoms under the coil, N = 2 p/(k T) pi r_w^2 l_coil (diatomic fill, 295 K). Exact.
- Minimum energy per ion pair: D0/2 + I + (3/2) Te (dissociation share, atomic ionization, electron heat at the
  PI temperature). Thermochemical constants; no radiation, no excitation, no wall loss. It is a strict lower
  bound on the real cost per ion pair E_c. A realistic E_c needs the H2/D2 reaction network (frc-transport.6),
  which is not built; E_c values above the minimum are shown as declared scan values, not computed ones.
- Upper bound on the ionized fraction: alpha <= eta W_PI / (N E_c), with eta the fraction of the stored PI energy
  delivered to the plasma (eta = 1 gives the hard bound).
- Flux freezing: the bias must stay frozen through the main rise. As in plasma_bench/frc_scaling.pi_temperature,
  the column skin time mu0 sigma r_w^2 must reach 10 main quarter periods. Here sigma = n_e e^2/(m_e (nu_ei+nu_en))
  includes electron-H2 momentum transfer: nu_en = n_n sigma_m <v_e>, sigma_m = 1.5e-19 m^2 (declared, with the
  doubled value as a sensitivity), and n_n = (1-alpha) N_atoms/2 undissociated molecules.
- PI bank scaling: the design window's ringing theta-PI table (150 kHz, 40 nH stray, 9.8 uF on the 77 nH coil)
  gives 60 J at B_PI = 0.028 T; stored energy, voltage and current scale as B_PI^2, B_PI and B_PI.
- Ionizing the remainder during the main implosion: the e-folding time of electron density in the H2 fill,
  1/(n_H2 K_iz(Te)), with K_iz from the repository's NIST-BEB total ionization of ground-state H2
  (plasma_bench/molecular_beb.py; Maxwellian, stationary targets). Atomic H and D2 are not computed (BEB here is
  H2-only; D2 is not aliased). Compared with the snowplow time, the main quarter period / 3 as in frc_scaling.
- Benchmark anchor: NUCTE-III (docs/NUCTE_BENCHMARK.md): D2 10 mTorr in a 0.23 m tube, 1.35 m coil,
  11.25 uF ringing at 124 kHz up to 135 kA, bias up to 0.09 T. Its stored PI energy follows from
  C, f and I alone (L = 1/(omega^2 C), V = I sqrt(L/C)); its field is the long-solenoid mu0 I/l.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.constants import mu_0, k as KB, e as QE, m_e, torr as TORR
from plasma_bench.molecular_beb import MolecularBEB

T_GAS = 295.
SIGMA_M = 1.5e-19                        # e-H2 momentum transfer near 1-10 eV [m^2], declared
GASES = dict(H2=dict(d0=4.478, ion=13.598), D2=dict(d0=4.556, ion=13.602))
DESIGN = dict(wall=.1524, coil_length=1.2, quarter_period=3.05e-6, bias_t=.05, bias_bank_j=5.3e3, main_bank_j=3.2e3)
PI_REF = dict(b=.028, energy=60., voltage=3.5e3, current=31e3, frequency=150e3)
FILLS = {'option A, H2 3.4 mTorr': ('H2', 3.4), 'option B, D2 3.5 mTorr': ('D2', 3.5)}


def atoms(fill_mtorr, radius, length):
    return 2*fill_mtorr*1e-3*TORR/(KB*T_GAS)*np.pi*radius**2*length


def minimum_cost_ev(gas, te_ev=2.):
    g = GASES[gas]
    return g['d0']/2+g['ion']+1.5*te_ev


def nucte_anchor():
    c, f, i = 11.25e-6, 124e3, 135e3
    inductance = 1/((2*np.pi*f)**2*c)
    voltage = i*np.sqrt(inductance/c)
    energy = .5*c*voltage**2
    n = atoms(10., .115, 1.35)
    return dict(atoms=n, energy_j=energy, voltage_v=voltage, inductance_nh=inductance*1e9,
                ev_per_atom=energy/(n*QE), b_pi_t=mu_0*i/1.35, bias_max_t=.09)


def skin_time(alpha, te_ev, n_atoms_m3, radius, sigma_m=SIGMA_M, coulomb_log=10.):
    n_e = alpha*n_atoms_m3
    nu_ei = 2.91e-12*n_e*coulomb_log*te_ev**-1.5
    nu_en = (1-alpha)*n_atoms_m3/2*sigma_m*np.sqrt(8*QE*te_ev/(np.pi*m_e))
    sigma = n_e*QE**2/(m_e*(nu_ei+nu_en))
    return mu_0*sigma*radius**2


def alpha_min(te_ev, n_atoms_m3, radius, required, sigma_m=SIGMA_M):
    """Smallest ionized fraction whose skin time reaches `required` (bisection in log alpha); None if unreachable."""
    if skin_time(1., te_ev, n_atoms_m3, radius, sigma_m) < required:
        return None
    lo, hi = -8., 0.
    for _ in range(80):
        mid = .5*(lo+hi)
        if skin_time(10**mid, te_ev, n_atoms_m3, radius, sigma_m) >= required:
            hi = mid
        else:
            lo = mid
    return 10**hi


def pi_bank(b_pi):
    s = b_pi/PI_REF['b']
    return dict(b_pi_t=b_pi, energy_j=PI_REF['energy']*s*s, voltage_v=PI_REF['voltage']*s, current_a=PI_REF['current']*s,
                below_bias=bool(b_pi <= DESIGN['bias_t']), bias_needed_t=max(b_pi, DESIGN['bias_t']),
                bias_bank_j=DESIGN['bias_bank_j']*(max(b_pi, DESIGN['bias_t'])/DESIGN['bias_t'])**2)


def budget():
    d = DESIGN
    anchor = nucte_anchor()
    out = dict(design=d, pi_reference=PI_REF, sigma_m_m2=SIGMA_M, nucte_iii=anchor, fills={})
    for name, (gas, mtorr) in FILLS.items():
        n = atoms(mtorr, d['wall'], d['coil_length'])
        n_m3 = 2*mtorr*1e-3*TORR/(KB*T_GAS)
        e_min = minimum_cost_ev(gas)
        row = dict(gas=gas, fill_mtorr=mtorr, atoms=n, atoms_m3=n_m3, minimum_cost_ev=e_min,
                   pi_ev_per_atom=PI_REF['energy']/(n*QE),
                   energy_full_ionization_at_minimum_j=n*e_min*QE,
                   alpha_bound={f'eta={eta:g}, Ec={ec:g} eV': min(1., eta*PI_REF['energy']/(n*ec*QE))
                                for eta in (1., .5, .25) for ec in (e_min, 30., 60., 120.)},
                   main_bank_fraction_to_ionize_rest={f'Ec={ec:g} eV': n*ec*QE/d['main_bank_j'] for ec in (e_min, 30., 60., 120.)},
                   flux_freezing={})
        required = 10*d['quarter_period']
        for te in (1., 2., 3., 5.):
            row['flux_freezing'][f'Te={te:g} eV'] = dict(
                alpha_min=alpha_min(te, n_m3, d['wall'], required),
                alpha_min_double_sigma_m=alpha_min(te, n_m3, d['wall'], required, 2*SIGMA_M),
                skin_time_full_us=skin_time(1., te, n_m3, d['wall'])*1e6)
        options = {}
        for label, target in (('dossier PI (as specified)', PI_REF['energy']/(n*QE)),
                              ('ideal full ionization at eta = 1', e_min),
                              ('NUCTE-III stored energy per atom', anchor['ev_per_atom'])):
            energy = target*n*QE
            b = PI_REF['b']*np.sqrt(energy/PI_REF['energy'])
            options[label] = dict(ev_per_atom=target, **pi_bank(b))
        options['dossier PI with the fill lowered to NUCTE-III eV/atom'] = dict(
            fill_mtorr=mtorr*PI_REF['energy']/(n*QE)/anchor['ev_per_atom'])
        row['pi_options'] = options
        if gas == 'H2':
            beb = MolecularBEB.from_nist('data/molecular/nist_beb', 'H2')
            snowplow = d['quarter_period']/3
            row['implosion_ionization'] = {f'Te={te:g} eV': dict(
                k_iz_m3_s=(k := beb.maxwellian_integral(te)['rate_within_energy_range_m3_s']),
                efold_us=1e6/(n_m3/2*k), efold_over_snowplow=1/(n_m3/2*k)/snowplow) for te in (5., 10., 20., 50.)}
        out['fills'][name] = row
    return out


def main(out='outputs/pi_budget'):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    result = budget()
    (out/'pi_budget.json').write_text(json.dumps(result, indent=1))
    a = result['nucte_iii']
    print(f"NUCTE-III anchor: {a['atoms']:.3g} atoms, PI {a['energy_j']:.0f} J ({a['voltage_v']/1e3:.1f} kV, "
          f"{a['inductance_nh']:.0f} nH), {a['ev_per_atom']:.0f} eV/atom, B_PI {a['b_pi_t']:.3f} T vs bias <= {a['bias_max_t']} T")
    for name, r in result['fills'].items():
        print(f"\n{name}: {r['atoms']:.3g} atoms under the coil; dossier PI {r['pi_ev_per_atom']:.1f} eV/atom; "
              f"minimum cost {r['minimum_cost_ev']:.2f} eV/ion; full ionization at the minimum {r['energy_full_ionization_at_minimum_j']:.0f} J")
        for k, v in r['alpha_bound'].items():
            print(f"  ionized fraction bound {k}: {v:.3f}")
        for k, v in r['flux_freezing'].items():
            print(f"  flux freezing {k}: alpha_min {v['alpha_min']} (2x sigma_m: {v['alpha_min_double_sigma_m']}), "
                  f"full-ionization skin time {v['skin_time_full_us']:.1f} us vs required {10*result['design']['quarter_period']*1e6:.1f} us")
        for k, v in r['main_bank_fraction_to_ionize_rest'].items():
            print(f"  main bank share to ionize the fill, {k}: {v:.3f}")
        for k, v in r.get('implosion_ionization', {}).items():
            print(f"  implosion ionization {k}: K_iz {v['k_iz_m3_s']:.2e} m3/s, e-fold {v['efold_us']:.2f} us "
                  f"= {v['efold_over_snowplow']:.2f} snowplow times")
        for k, v in r['pi_options'].items():
            print(f"  option '{k}': {json.dumps({q: (round(x, 4) if isinstance(x, float) else x) for q, x in v.items()})}")
    return result


if __name__ == '__main__':
    main(*sys.argv[1:2])
