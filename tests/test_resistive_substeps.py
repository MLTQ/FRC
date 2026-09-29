"""Stiff-decay control against an independently diagonalized exact solution."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from scipy.linalg import eigh
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.material_resistivity import MaterialResistivity
from plasma_bench.resistive_substeps import integrate_diffusion


def fixture():
    d=MovingDomain(4,8);o=MaterialResistivity(d,d.fluid_reference,.1*MU0)
    n=d.fluid_node_count;k=o.stiffness.toarray();m=o.mass.toarray()
    s=k[:n,:n]-k[:n,n:]@np.linalg.solve(k[n:,n:],k[n:,:n])
    values,vectors=eigh(s,m)
    initial=vectors[:,0]*np.sqrt(MU0)/vectors[0,0]
    return o,initial,4/(.1*values[0]),m


def test_adaptive_diffusion_resolves_raw_midpoint_false_reversal():
    o,a,duration,m=fixture();raw=o.step(a,duration)
    assert (a@m@raw.material_a)/(a@m@a)<0
    errors=[]
    for tolerance in (1e-3,1e-4,1e-5):
        result=integrate_diffusion(o,a,duration,tolerance=tolerance)
        amplitude=(a@m@result.material_a)/(a@m@a)
        errors.append(abs(amplitude-np.exp(-4)))
        assert amplitude>0
        assert result.rejected_trials>0 and result.accepted_substeps>2
        assert result.maximum_accepted_error<=tolerance
        final=o.induced(result.material_a)
        energy=np.pi/MU0*(final@(o.stiffness@final))
        assert abs(energy+result.cell_heat.sum()-raw.initial_energy)<1e-11
    assert errors[-1]<2e-5
    assert np.all(np.array(errors[:-1])/errors[1:] > 2.5)


def test_diffusion_trial_exhaustion_leaves_input_unchanged():
    o,a,duration,_=fixture();before=a.copy()
    small=integrate_diffusion(o,a,1e-6,tolerance=1e-4,max_trials=1)
    assert small.accepted_substeps==2 and small.rejected_trials==0
    with pytest.raises(RuntimeError,match="trial budget"):
        integrate_diffusion(o,a,duration,tolerance=1e-8,max_trials=1)
    assert np.array_equal(a,before)
