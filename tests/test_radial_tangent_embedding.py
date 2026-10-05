"""Independent one-sided analytic polynomial subdivision."""
from types import SimpleNamespace
import numpy as np
from examples.radial_tangent_embedding import embed


def grid(elements):
    radius=2.;degree=5;nz=2
    # Arbitrary distinct nodes suffice: the interpolation contract is algebraic.
    edges=.5*np.linspace(0,radius,elements+1)**2
    x=(edges[:-1,None]+np.diff(edges)[:,None]*np.linspace(0,1,degree+1))[...,None]
    return SimpleNamespace(elements=elements,nz=nz,degree=degree,radius=radius,x=x,shape=(elements,degree+1,nz,6))


def test_exact_discontinuous_polynomial_subdivision():
    old=grid(2);new=grid(8);value=np.zeros(old.shape);expected=np.zeros(new.shape)
    for mesh,array in ((old,value),(new,expected)):
        for child in range(mesh.elements):
            parent=child*old.elements//mesh.elements;x=mesh.x[child,:,0]
            for k in range(6):
                power=4 if parent==0 and k in (1,2) else 5
                array[child,:,:,k]=(7*parent+k+3*x**power)[:,None]*np.array([1.,2.])
        array[0,0,:,1:3]=0.
    actual=embed(old,new,value)
    np.testing.assert_allclose(actual,expected,rtol=3e-14,atol=3e-13)
    # The original interface jump must remain one-sided after subdivision.
    np.testing.assert_allclose(actual[4,0]-actual[3,-1],expected[4,0]-expected[3,-1],atol=3e-13)
