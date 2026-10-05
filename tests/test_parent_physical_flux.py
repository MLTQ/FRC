"""Independent expanded finite-mass kinetic flux and Poynting controls."""
import numpy as np
from scipy.constants import elementary_charge as E,electron_mass as ME,atomic_mass,mu_0 as MU0
from examples.physical_parent_flux import fluxes


class UniformSampler:
    def __init__(self,n,v,entropy):
        self.n=np.array([[n]]);self.v=np.asarray(v)[None,None,:];self.entropy=np.asarray(entropy)
        self.d=type('D',(),{'thermal':type('T',(),{'n0':n,'gammas':np.array([5/3,1.4])})()})()
    def fluid(self,state):return self.n,self.v,None
    def scalar(self,value):return np.array([[value]])
    def electric(self,value):return np.asarray(value)[None,None,:]
    def magnetic(self,value):return np.asarray(value)[None,None,:]


def test_parent_flux_matches_independent_expanded_mass_moments():
    n=3e19;v=np.array([250.,-310.,750.]);j=np.array([170.,250.,-40.]);electric=np.array([4.,-8.,2.])
    magnetic=np.array([.003,.002,-.001]);bias=.006;mi=39.948*atomic_mass
    entropy=np.log([3.,.7]);state=np.r_[n,np.zeros(3),n*entropy]
    s=UniformSampler(n,v,entropy);f=fluxes(s,state,magnetic,electric,j,bias,mi)[0,0]
    mass=mi+ME;mu=mi*ME/mass;delta=(mi-ME)/mass;w=j/(E*n)
    kinetic=.5*mass*n*np.dot(v,v)*v+.5*mu*n*np.dot(w,w)*v+mu*n*np.dot(v,w)*w-.5*mu*delta*n*np.dot(w,w)*w
    gamma=s.d.thermal.gammas;h=E*np.exp(entropy)/(gamma-1);a=np.array([-mi/mass,ME/mass])
    enthalpy=np.sum(gamma*n*h)*v+np.sum(a*gamma*h)/E*j
    b=magnetic.copy();b[2]+=bias
    expected=kinetic+enthalpy+np.cross(electric,b)/MU0
    assert np.allclose(f[-1],expected,rtol=3e-15,atol=1e-12)
    assert np.allclose(f[0],n*v,rtol=1e-15)
    assert np.array_equal(f[1],j)
    assert np.allclose(f[2],n*entropy[0]*v-mi/mass*entropy[0]/E*j,rtol=3e-15)
