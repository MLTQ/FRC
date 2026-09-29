"""Independent spatial eigenmode audit of current-inertia time integration."""
import numpy as np
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from scipy.linalg import eigh
from plasma_bench.current_inertia import MaterialCurrent
from plasma_bench.moving_geometry import MovingDomain


def measure():
    d=MovingDomain(4,8);de=.2;diffusivity=.1
    o=MaterialCurrent(d,d.fluid_reference,electron_density=ME/(MU0*E**2*de**2),resistivity=diffusivity*MU0)
    n=d.fluid_node_count;k=o.stiffness.toarray();m=o.mass.toarray()
    schur=k[:n,:n]-k[:n,n:]@np.linalg.solve(k[n:,n:],k[n:,:n]);values,vectors=eigh(schur,m)
    mode=5;initial=vectors[:,mode]*1e-8;rate=diffusivity*values[mode]/(1+de*de*values[mode]);end=.2
    exact=initial*np.exp(-rate*end);rows=[]
    for dt in (.04,.02,.01):
        a=initial.copy()
        for _ in range(round(end/dt)):step=o.step(a,dt);a=step.material_a
        error=float(np.sqrt((a-exact)@m@(a-exact)/(exact@m@exact)))
        rows.append({'dt_s':dt,'relative_current_L2_error':error,'final_step_balance_residual_j':step.balance_residual})
    return {'configuration':{'radial_cells':4,'angular_cells':8,'radius_m':1.,'matching_radius_m':2.,
              'modes':8,'quadrature':5,'electron_length_m':de,'diffusivity_m2_s':diffusivity,'end_time_s':end,
              'eigenmode_zero_based':mode,'initial_amplitude_multiplier_t':1e-8},
            'interpretation':'Independent dense spatial eigensolve isolates time order; not a continuum spatial reference',
            'spatial_eigenvalue_m2_inverse':float(values[mode]),'exact_decay_rate_s_inverse':float(rate),
            'physical_high_wavenumber_rate_limit_s_inverse':diffusivity/de**2,'cases':rows,
            'error_halving_ratios':[rows[i]['relative_current_L2_error']/rows[i+1]['relative_current_L2_error'] for i in range(2)]}
