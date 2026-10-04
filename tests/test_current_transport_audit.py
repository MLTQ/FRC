"""Analytic checks of the independent physical transport audit."""
import numpy as np
from examples.current_transport_sampling import TransportSampler,bases
from examples.regular_coupling_reference import initialize


def test_root_edge_derivatives_match_exact_polynomial_and_periodic_boundary():
    edges=np.linspace(-1,1,9);x=np.array([-.97,-.2,.13,.9]);n,d,dd,e,de=bases(edges,x,5)
    assert np.allclose(n@edges**4,x**4,rtol=1e-12,atol=1e-13)
    assert np.allclose(d@edges**4,4*x**3,rtol=1e-12,atol=1e-13)
    assert np.allclose(dd@edges**4,12*x*x,rtol=1e-12,atol=1e-13)
    average=np.diff(edges**4)/(4*np.diff(edges))
    assert np.allclose(e@average,x**3,rtol=1e-12,atol=1e-13)
    assert np.allclose(de@average,3*x*x,rtol=1e-12,atol=1e-13)
    left=bases(edges,[-1+1e-7],5,True);right=bases(edges,[1+1e-7],5,True)
    for a,b in zip(left,right):assert np.allclose(a,b,rtol=2e-12,atol=2e-13)


def test_independent_velocity_rate_is_derivative_of_declared_physical_fields():
    ref,d,u,m,_=initialize(dict(elements=2));s=TransportSampler(d,m.mesh,8);rate=np.zeros_like(u)
    rate[...,0]=u[...,0]*(.2+.1*(d.r/ref.radius)**2)
    rate[...,1]=.3*u[...,1];rate[...,2]=-.1*u[...,2];rate[...,3]=.4*u[...,3]
    n,v,grad,nt,vt=s.fluid(u,rate);h=1e-3
    # Richardson cancellation checks the quotient derivative without reusing it.
    def difference(step):return (s.fluid(u+step*rate)[1]-s.fluid(u-step*rate)[1])/(2*step)
    numerical=(4*difference(h/2)-difference(h))/3
    assert np.linalg.norm(numerical-vt)/np.linalg.norm(vt)<3e-11
    assert np.allclose(nt,s.scalar(rate[...,0]))
