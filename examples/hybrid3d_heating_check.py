"""Complete resistive hybrid runs against force-free magnetic decay and electron heating."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE, m_p as MP
from scipy.integrate import solve_ivp
from plasma_bench.hybrid3d import Box, HybridParams
from plasma_bench.hybrid3d_device import Device, DeviceHybrid


def decay(cells, resistivity, backend='cpu', gpu_moments=False, fused_arrays=False):
    """B=A(cos kz,sin kz,0) is force-free: ions stay still, heat stays uniform.

    A'=-eta(Te) k² A/mu0 and 3n e (Te-Te0)=(A0²-A²)/mu0.
    The reference uses continuum k; no production operator enters its ODE.
    """
    length, amplitude, density, te0 = .2, .1, 1e22, 5.
    k, gamma = 2*np.pi/length, 5/3
    eta0 = 1e-3 if resistivity == 'constant' else 1.03e-3*te0**-1.5
    duration = .5*MU0/(eta0*k*k)
    dt = min(.05*MP/(QE*amplitude), duration/100)
    steps = int(np.ceil(duration/dt)); dt = duration/steps
    dz = length/cells
    box = Box((2, 2, cells), (0., 0., 0.), (2*dz, 2*dz, length))
    phase = k*(np.arange(cells)+.5)*dz
    b = [np.broadcast_to(amplitude*np.cos(phase), box.shape).copy(),
         np.broadcast_to(amplitude*np.sin(phase), box.shape).copy(), np.zeros(box.shape)]
    points = np.meshgrid(*[(np.arange(q)+.5)*box.d[a] for a, q in enumerate(box.shape)], indexing='ij')
    x = np.column_stack([p.ravel() for p in points])
    model = DeviceHybrid(box, Device(np.inf, np.inf, np.inf, 0.),
                         HybridParams(ion_mass=MP, te_ev=te0, gamma_e=gamma, n_floor=.01*density),
                         b, x, np.zeros_like(x), density*box.cell_volume, dt, n_ref=density,
                         thermal=dict(gamma=gamma, resistivity=resistivity, constant_eta=eta0,
                                      backend=backend, gpu_moments=gpu_moments,
                                      fused_arrays=fused_arrays))
    initial_energy = model.magnetic_energy()+model.electrons.energy(model.n)
    for _ in range(steps):
        model.step()
    measured = float(np.mean(model.b[0]*np.cos(phase)+model.b[1]*np.sin(phase)))

    def temperature(a):
        return te0+(amplitude**2-a*a)/(3*MU0*density*QE)

    def rhs(t, y):
        eta = eta0 if resistivity == 'constant' else 1.03e-3*max(temperature(y[0]), 1.)**-1.5
        return [-eta*k*k*y[0]/MU0]
    reference = solve_ivp(rhs, (0., duration), [amplitude], method='DOP853', rtol=1e-12, atol=1e-14)
    exact = float(reference.y[0, -1])
    te = float(np.mean(model.electrons.temperature(model.n)))
    total = model.magnetic_energy()+model.electrons.energy(model.n)+model.kinetic_energy()
    wave_energy = amplitude**2*np.prod(np.asarray(box.hi)-box.lo)/(2*MU0)
    return dict(cells=cells, resistivity=resistivity, duration_s=duration, steps=steps, substeps=model.substeps,
                B_exact_T=exact, B_measured_T=measured, B_relative_error=(measured-exact)/exact,
                Te_exact_eV=temperature(exact), Te_measured_eV=te,
                total_energy_error_over_initial_magnetic=(total-initial_energy)/wave_energy,
                ion_kinetic_j=model.kinetic_energy(), cutoff_removed_heat_j=model.electrons.cutoff_removed_heat_j,
                div_B_relative=model.max_divergence())


def main(out='outputs/hybrid3d_heating_check'):
    out = Path(out); out.mkdir(exist_ok=True, parents=True)
    rows = []
    for resistivity in ('constant', 'spitzer'):
        for cells in (16, 32, 64):
            row = decay(cells, resistivity); rows.append(row)
            print(json.dumps(row), flush=True)
            (out/'comparison.json').write_text(json.dumps(rows, indent=1))
    return rows


def entropy_contact(cells, backend='cpu', gpu_moments=False, fused_arrays=False):
    """Uniform pressure, sinusoidal density, constant drift: exact translation for one period."""
    length, density, te0, velocity, amplitude = .2, 1e21, 5., 1e5, .3
    dx = length/cells
    box = Box((cells, 2, 2), (0., 0., 0.), (length, 2*dx, 2*dx))
    xyz = np.meshgrid(*[(np.arange(q)+.5)*dx for q in box.shape], indexing='ij')
    x = np.column_stack([p.ravel() for p in xyz])
    w = density*(1+amplitude*np.cos(2*np.pi*x[:, 0]/length))*box.cell_volume
    v = np.zeros_like(x);v[:, 0] = velocity
    # CIC node density of this lattice is the adjacent-centre average.
    xn = np.arange(cells)*dx
    initial_te = te0/(1+amplitude*np.cos(2*np.pi*xn/length)*np.cos(np.pi/cells))[:, None, None]
    steps = 5*cells;dt = length/(velocity*steps)
    model = DeviceHybrid(box, Device(np.inf, np.inf, np.inf, 0.),
                         HybridParams(ion_mass=MP, te_ev=te0, gamma_e=5/3, n_floor=.01*density),
                         [np.zeros(box.shape) for _ in range(3)], x, v, w, dt, n_ref=density,
                         thermal=dict(resistivity='constant', constant_eta=0., initial_temperature_ev=initial_te,
                                      backend=backend, gpu_moments=gpu_moments,
                                      fused_arrays=fused_arrays))
    energy0 = model.electrons.energy(model.n)+model.kinetic_energy()
    heat0 = model.electrons.energy(model.n)
    for _ in range(steps):model.step()
    pressure = model._pressure(model.n)
    exact = density*te0*QE
    return dict(cells=cells, pressure_L1_relative=float(np.mean(abs(pressure-exact))/exact),
                max_velocity_relative_error=float(max(np.max(abs(model.v[:, 0]-velocity)),
                                                      np.max(abs(model.v[:, 1:])))/velocity),
                total_energy_error_over_initial_electron=(model.electrons.energy(model.n)+model.kinetic_energy()-energy0)/heat0)


def energy_advection(cells, backend='cpu', fused_arrays=False):
    """Electron transport operator at prescribed n,u; exact periodic translation.

    This deliberately isolates the PDE operator, not a force-balanced device.
    An outer CFL near1.5 requires the new positive-energy substeps.
    """
    length, density, speed = .2, 1e20, 1e5
    dx = length/cells
    box = Box((cells,2,2),(0.,0.,0.),(length,2*dx,2*dx))
    xyz = np.meshgrid(*[(np.arange(q)+.5)*dx for q in box.shape],indexing='ij')
    x = np.column_stack([p.ravel() for p in xyz])
    te = (1e-8+5*np.exp(-.5*((np.arange(cells)*dx-.5*length)/(.04*length))**2))[:,None,None]
    model = DeviceHybrid(box,Device(np.inf,np.inf,np.inf,0.),
        HybridParams(ion_mass=MP,te_ev=5.,gamma_e=5/3,n_floor=.01*density),
        [np.zeros(box.shape) for _ in range(3)],x,np.zeros_like(x),density*box.cell_volume,
        1e-9,n_ref=density,thermal=dict(backend=backend,resistivity='constant',constant_eta=0.,
            initial_temperature_ev=te,fused_arrays=fused_arrays))
    electrons=model.electrons
    electrons.bulk_velocity[0].fill(speed)
    reference=electrons.u.copy();energy0=electrons.energy(model.n)
    steps=int(np.ceil(cells/1.5));dt=length/(speed*steps)
    zeros=np.zeros_like(model.n);vectors=np.zeros((3,)+model.n.shape)
    for _ in range(steps):
        electrons.advance(model.n,vectors,vectors,dt,zeros,hall_work=zeros,hall_velocity=vectors)
        model.time+=dt
    return dict(cells=cells,backend=backend,outer_cfl=speed*dt/dx,
        relative_L1=float(np.mean(abs(electrons.u-reference))/np.mean(reference)),
        relative_energy_error=(electrons.energy(model.n)-energy0)/energy0,
        minimum_U=float(electrons.u.min()),max_subcycles=electrons.max_energy_subcycles,
        actual_max_loss_fraction=electrons.max_energy_loss_fraction,
        roundoff_energy_j=electrons.roundoff_energy_j)


if __name__ == '__main__':
    main()
