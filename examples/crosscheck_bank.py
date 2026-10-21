"""Independent three-branch lumped main-bank calculation; no project circuit imports."""
import json
import sys
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class BankCase:
    name: str
    coil_h: float
    target_a: float
    bias_a: float
    dossier_kv: float
    dossier_peak_ka: float
    dossier_slew_ka_us: float
    dossier_crowbar_us: float
    dossier_reverse_kv: float


CASES = (
    BankCase('A, 0.2 T', 77e-9, 221e3, -55.25e3, 11.6, 36, 18, 3.05, 2.9),
    BankCase('A, 0.3 T', 77e-9, 332e3, -55.25e3, 16.2, 50, 26, 3.05, 4.0),
    BankCase('B, 0.2 T', 93e-9, 266e3, -66.5e3, 15.2, 43, 20, 3.33, 3.5),
    BankCase('B, 0.3 T', 93e-9, 399e3, -66.5e3, 21.3, 60, 29, 3.33, 4.8),
)

N, C, LS, LCB, RCB = 8, 48e-6, 30e-9/8, 8e-9/8, 2.1e-3/8
CBIAS, LBIAS, RBIAS = 2e-3, 2.4e-6, 8e-3


def integrate(case, voltage, crowbar=False, initial=None):
    """One switch interval: L dI/dt = V-RI, dV/dt = -I/C.

    Coil inductance links every active branch; diagonal inductances are the
    bias isolation, parallel main leads and parallel crowbar loops. Branch
    currents add at the coil, with the bias negative. No plasma load is used.
    """
    active = 3 if crowbar else 2
    inductance = np.diag([LBIAS, LS, LCB][:active])+case.coil_h*np.ones((active, active))
    resistance = np.array([RBIAS, 0., RCB][:active])
    inverse = np.linalg.inv(inductance)
    y0 = np.zeros(active+3)
    if initial is None:
        y0[0] = case.bias_a
        y0[active:active+2] = [RBIAS*case.bias_a, voltage]
    else:
        y0[:2] = initial[:2]
        y0[active:] = initial[2:]

    def rhs(t, y):
        v = np.zeros(active); v[:2] = y[active:active+2]
        current = y[:active]
        return np.r_[inverse@(v-resistance*current), -current[0]/CBIAS,
                     -current[1]/C, np.sum(resistance*current**2)]

    def event(t, y):
        return y[1] if crowbar else y[active+1]
    event.direction, event.terminal = -1, True
    result = solve_ivp(rhs, (0., 10e-6), y0, method='DOP853', rtol=2e-11,
                       atol=1e-8, max_step=1e-8, events=event, dense_output=True)
    if not result.success or not len(result.t_events[0]):
        raise RuntimeError('Switch event was not reached: '+result.message)
    return result, inductance, resistance, rhs(0., y0)


def calculate(case, voltage):
    before, l0, _, slope = integrate(case, voltage)
    after, l1, _, _ = integrate(case, voltage, True, before.y[:, -1])
    times0 = np.linspace(0., before.t[-1], 2001)
    times1 = np.linspace(0., after.t[-1], 1001)
    a, b = before.sol(times0), after.sol(times1)
    peak = max(float(a[1].max()), float(b[1].max()))
    coil = max(float(a[:2].sum(axis=0).max()), float(b[:3].sum(axis=0).max()))

    def energy(y, inductance):
        n = inductance.shape[0]
        return (.5*np.einsum('it,ij,jt->t', y[:n], inductance, y[:n])
                +.5*CBIAS*y[n]**2+.5*C*y[n+1]**2+y[n+2])
    initial = energy(before.y[:, :1], l0)[0]
    residual = max(float(abs(energy(a, l0)-initial).max()), float(abs(energy(b, l1)-initial).max()))
    quarter = .5*np.pi*np.sqrt(C*(LS+case.coil_h*LBIAS/(case.coil_h+LBIAS)))
    return dict(charge_kv=voltage/1e3, bank_energy_kj=.5*C*voltage**2/1e3,
                coil_peak_ka=coil/1e3, switch_peak_ka=peak/(N*1e3),
                initial_slew_ka_us=slope[1]/(N*1e9), crowbar_us=before.t[-1]*1e6,
                capacitor_reverse_kv=-float(after.y[4, -1])/1e3,
                switch_open_us=(before.t[-1]+after.t[-1])*1e6,
                analytic_quarter_us=quarter*1e6, energy_residual_j=residual)


def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in CASES:
        v = brentq(lambda voltage: calculate(case, voltage)['coil_peak_ka']*1e3-case.target_a, 5e3, 30e3, xtol=1e-6)
        actual = calculate(case, v)
        reference = dict(charge_kv=case.dossier_kv, switch_peak_ka=case.dossier_peak_ka,
                         initial_slew_ka_us=case.dossier_slew_ka_us, crowbar_us=case.dossier_crowbar_us,
                         capacitor_reverse_kv=case.dossier_reverse_kv)
        row = dict(case=case.name, independent=actual, dossier=reference,
                   difference_percent={k:100*(actual[k]/value-1) for k,value in reference.items()})
        rows.append(row); print(json.dumps(row), flush=True)
    (out/'main_bank.json').write_text(json.dumps(rows, indent=2)+'\n')
    return rows


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv)>1 else 'outputs/crosscheck_bank')
