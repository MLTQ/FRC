"""Independent analytic cylindrical power and Fourier localization."""
from types import SimpleNamespace
import numpy as np
from examples.extracted_response_fields import localization


def test_cylindrical_polynomial_fourier_power():
    q,w=np.polynomial.legendre.leggauss(20);nr=8;nz=16;radius=.3;length=.8
    rf=np.linspace(0,radius,nr+1);zf=np.linspace(-length/2,length/2,nz+1)
    r=(.5*(rf[:-1]+rf[1:])[:,None]+radius/(2*nr)*q).ravel()
    z=(.5*(zf[:-1]+zf[1:])[:,None]+length/(2*nz)*q).ravel()
    wr=(np.pi*radius/nr*r.reshape(nr,-1)*w).ravel()
    wz=np.broadcast_to(length/(2*nz)*w,(nz,len(q))).ravel()
    sample=SimpleNamespace(r=r,z=z,weight=wr[:,None]*wz[None,:],
        mesh=SimpleNamespace(nz=nz,rf=rf),d=SimpleNamespace(radius=radius,length=length,elements=4))
    value=r[:,None,None]*np.stack((np.ones_like(z),np.cos(6*np.pi*z/length),np.zeros_like(z)),axis=-1)[None,:,:]
    result,raw=localization(sample,value)
    # Integral 2*pi*r^3 dr times L*(1+1/2), with exact Fourier k=0,+/-3.
    total=3*np.pi*length*radius**4/4
    np.testing.assert_allclose(raw['fourier_power'].sum(),total,rtol=2e-14)
    np.testing.assert_allclose(result['component_power_fractions'],[2/3,1/3,0],atol=2e-14)
    np.testing.assert_allclose(result['first_element_power_fraction'],1/4**4,atol=2e-14)
    np.testing.assert_allclose(result['last_element_power_fraction'],1-(3/4)**4,atol=2e-14)
    np.testing.assert_allclose(result['outer_field_cell_power_fraction'],1-(7/8)**4,atol=2e-14)
    expected=np.zeros_like(raw['fourier_power']);k=raw['fourier_index']
    expected[k==0]=2*total/3;expected[abs(k)==3]=total/6
    np.testing.assert_allclose(raw['fourier_power'],expected,atol=1e-15)
