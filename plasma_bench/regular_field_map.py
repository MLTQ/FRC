"""Physical samples and force-work adjoints of the actual regular fluid state."""
import numpy as np
from scipy.sparse import csr_matrix


def polynomial(nodes, points):
    """Lagrange products; only small local element matrices are formed."""
    result=np.ones((len(points),len(nodes)))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if i!=j:result[:,i]*=(points-nodes[j])/(nodes[i]-nodes[j])
    return result


def periodic(nodes, points, length):
    """Real periodic interpolant, including the cosine Nyquist mode when even."""
    theta=2*np.pi*(np.asarray(points)[:,None]-np.asarray(nodes)[None,:])/length
    n=len(nodes);result=np.ones_like(theta)
    for k in range(1,(n-1)//2+1):result+=2*np.cos(k*theta)
    if n%2==0:result+=np.cos((n//2)*theta)
    return result/n


class RegularFieldMap:
    def __init__(self,fluid,geometry):
        mesh=geometry.mesh
        if (not np.isclose(mesh.radius,fluid.radius,rtol=1e-14,atol=0.) or
            not np.isclose(2*mesh.half_length,fluid.length,rtol=1e-14,atol=0.) or
            mesh.nr%fluid.elements):
            raise ValueError('Matching cylinder and field cells subdividing fluid elements required')
        self.fluid=fluid;self.geometry=geometry;self.shape=fluid.shape[:-1]
        r=geometry.radial_points;elements=np.minimum((r*fluid.elements/fluid.radius).astype(int),fluid.elements-1)
        maps=[]
        for transverse in (False,True):
            matrix=np.zeros((len(r),fluid.elements*(fluid.degree+1)))
            for e in range(fluid.elements):
                rows=np.flatnonzero(elements==e);begin=1 if transverse and e==0 else 0
                columns=e*(fluid.degree+1)+np.arange(begin,fluid.degree+1)
                matrix[np.ix_(rows,columns)]=polynomial(fluid.x[e,begin:,0],.5*r[rows]**2)
            maps.append(csr_matrix(matrix))
        self.radial,self.radial_velocity=maps
        self.axial=csr_matrix(periodic(fluid.z.ravel(),geometry.axial_points,fluid.length))

    def validate(self,values):
        values=np.asarray(values,float)
        if values.shape!=self.shape or not np.isfinite(values).all():raise ValueError('Finite regular scalar coefficients required')
        return values

    def scalar(self,values,*,transverse=False):
        radial=self.radial_velocity if transverse else self.radial
        result=self.geometry.tensor(radial,self.validate(values).reshape(radial.shape[1],self.fluid.nz),self.axial)
        if not np.isfinite(result).all():raise ValueError('Nonfinite physical regular samples')
        return result

    def adjoint(self,values,*,transverse=False):
        geo=self.geometry;values=np.asarray(values,float)
        if values.shape!=geo.r.shape or not np.isfinite(values).all():raise ValueError('Finite scalar at physical quadrature required')
        radial=self.radial_velocity if transverse else self.radial
        weighted=geo._unlayout(geo.weight*values)
        result=np.asarray((radial.T@weighted)@self.axial).reshape(self.shape)
        if not np.isfinite(result).all():raise ValueError('Nonfinite physical regular adjoint')
        return result

    def velocity(self,state):
        self.fluid.fields(state);n=self.scalar(state[...,0])
        if np.any(n<=0):raise ValueError('Positive reconstructed physical density required')
        return np.stack((self.geometry.r*self.scalar(state[...,1]/state[...,0],transverse=True),
                         self.geometry.r*self.scalar(state[...,2]/state[...,0],transverse=True),
                         self.scalar(state[...,3])/n),axis=-1)

    def force_rate(self,state,force):
        self.fluid.fields(state);geo=self.geometry;d=self.fluid;force=np.asarray(force,float)
        if force.shape!=geo.r.shape+(3,) or not np.isfinite(force).all():raise ValueError('Finite physical force density required')
        n=self.scalar(state[...,0])
        if np.any(n<=0):raise ValueError('Positive reconstructed physical density required')
        rate=np.zeros(d.shape);denominator=d.mass*d.volume*2*d.x;active=denominator>0
        for k in (0,1):
            weak=self.adjoint(geo.r*force[...,k],transverse=True)
            np.divide(weak,denominator,out=rate[...,k+1],where=active)
        rate[...,3]=state[...,0]*self.adjoint(force[...,2]/n)/(d.mass*d.volume)
        if not np.isfinite(rate).all():raise ValueError('Nonfinite regular force rate')
        return rate


class RegularDensity:
    """Immutable actual regular point density and its common physical pullback."""
    def __init__(self,mapping,number_density):
        self.mapping=mapping;self.geometry=mapping.geometry
        self.point_values=mapping.validate(number_density).copy();self.samples=self.sample(self.point_values)
        if np.any(self.point_values<=0) or np.any(self.samples<=0):raise ValueError('Positive point and sampled density required')
        self.point_values.setflags(write=False);self.samples.setflags(write=False)

    def sample(self,values):return self.mapping.scalar(values)

    def adjoint(self,values):return self.mapping.adjoint(values)
