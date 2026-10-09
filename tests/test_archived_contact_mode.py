"""Independent cylindrical curl of the fixed high-radial reference."""
import numpy as np
from scipy.constants import mu_0 as MU0
from examples.archived_contact_mode import ArchivedContactMode


def test_high_radial_reference_strong_ampere_by_finite_differences():
    mode=ArchivedContactMode();r=np.array([.009,.017,.029]);z=np.array([-.083,-.014,.071])
    def magnetic(rr,zz):return mode.fields(rr,zz)['physical'][...,6:9]
    def derivative(axis,h):
        def f(offset):return magnetic(r+offset*h if axis==0 else r,z+offset*h if axis==2 else z)
        return (f(-2)-8*f(-1)+8*f(1)-f(2))/(12*h)
    dr=derivative(0,1e-6);dz=derivative(2,1e-5);b=magnetic(r,z)
    curl=np.stack((-dz[...,1],dz[...,0]-dr[...,2],dr[...,1]+b[...,1]/r),axis=-1)
    expected=mode.fields(r,z)['current']
    for k in range(3):
        error=np.linalg.norm(curl[...,k]/MU0-expected[...,k])/np.linalg.norm(expected[...,k])
        assert error<2e-8
