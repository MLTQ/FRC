"""Analytic cylindrical induction and original strict/atomic controls."""
import numpy as np
import pytest
from scipy.constants import mu_0 as MU0
from examples.strong_wall_reference import initialize
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.variable_compatible_current import weak_vector
from plasma_bench.split_inverse_wall_induction import SplitInverseWallInduction


@pytest.mark.parametrize('nr',[6,12,24])
def test_analytic_polynomial_induction_full_fields_and_work(nr):
    _,_,u,s,_=initialize(dict(kind='contact',elements=2,points=8,field_nr=nr,field_nz=8))
    model=RegularCurrent(s.mapping,u[...,0],s.ion_mass);geo=s.geometry;r=geo.mesh.radius;x=geo.r/r
    electric=np.zeros(geo.r.shape+(3,));electric[...,1]=x*(2-3*x*x+x**4);electric[...,2]=.37
    current=np.zeros_like(electric);current[...,1]=24*x*(x*x-1)/(MU0*r*r)
    magnetic=np.zeros_like(electric);magnetic[...,2]=(-4+12*x*x-6*x**4)/r
    forcing=weak_vector(geo,electric-model.inertia_coefficient[...,None]*current)
    solver=SplitInverseWallInduction(model,tolerance=1e-10);result=solver.solve(forcing)
    actual=(geo.sample_electric(result.electric),geo.sample_electric(result.current),geo.magnetic_samples(*geo.unpack(result.magnetic)))
    for value,expected in zip(actual,(electric,current,magnetic)):
        assert np.linalg.norm(value-expected)/np.linalg.norm(expected)<2e-10
    assert result.residual<=1e-10
    # Original rate-work identity: -Jt.f equals Bt.H.Bt + Jt.I.Jt.
    powers=[-result.current@forcing,result.magnetic@(geo.hessian@result.magnetic),result.current@(model.inertia@result.current)]
    assert abs(powers[0]-powers[1]-powers[2])/sum(abs(v) for v in powers)<1e-10


def test_exhaustion_leaves_forcing_and_actual_density_unchanged():
    _,_,u,s,_=initialize(dict(kind='contact',elements=2,points=8,field_nr=6,field_nz=8))
    model=RegularCurrent(s.mapping,u[...,0],s.ion_mass)
    f=np.random.default_rng(156).normal(size=model.mass.shape[0]);before=f.copy();density=model.density.copy()
    with pytest.raises(RuntimeError,match='residual failure'):
        SplitInverseWallInduction(model,tolerance=1e-10,max_iterations=1).solve(f)
    assert np.array_equal(f,before) and np.array_equal(model.density,density)
