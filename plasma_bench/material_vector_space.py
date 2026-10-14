"""Full axisymmetric electric/current and solenoidal material magnetic spaces."""
import numpy as np
from scipy.linalg import null_space
from scipy.sparse import csr_matrix,hstack,kron
from .compatible_polynomials import bases
from .material_current_space import MaterialCurrentSpace


def trace_kernel(trace):
    """Condition rows of mixed value/derivative traces without changing their null space."""
    scale=np.linalg.norm(trace,axis=1)
    if np.any(scale==0):raise ValueError('Independent nonzero physical trace rows required')
    return csr_matrix(null_space(trace/scale[:,None]))


class MaterialVectorSpace:
    """Three physical components, including a regular poloidal flux coordinate.

    This first vector boundary contract has B=0 at the material surface.
    Penetrating vacuum/coil flux is not supplied by this compact-field space.
    """
    def __init__(self,mesh,*,degree=3,quadrature=6):
        base=MaterialCurrentSpace(mesh,degree=degree,quadrature=quadrature)
        self.mesh=mesh;self.degree=degree;self.order=quadrature
        self.reference_r=base.reference_r;self.reference_z=base.reference_z
        q=quadrature;r=base.reference_r.reshape(mesh.nr,mesh.nz,q,q)[:,0,:,0].ravel()
        z=base.reference_z.reshape(mesh.nr,mesh.nz,q,q)[0,:,0,:].ravel()
        nx,dx,_,_=bases(.5*mesh.rf**2,.5*r*r,degree)
        nz,dz,_,_=bases(mesh.zf,z,degree)
        order=np.arange(len(r)*len(z)).reshape(mesh.nr,q,mesh.nz,q).transpose(0,2,1,3).ravel()
        scalar=kron(nx,nz,format='csr')[order]
        toroidal=kron(nx.multiply(r[:,None]),nz,format='csr')[order]
        curl_r=-kron(nx.multiply(r[:,None]),dz,format='csr')[order]
        curl_z=kron(2*nx+dx.multiply((r*r)[:,None]),nz,format='csr')[order]
        points=len(order);self.scalar_count=scalar.shape[1];self.meridional_electric_count=base.electric_count
        zero_scalar=csr_matrix((points,self.scalar_count));zero_meridional=csr_matrix((points,base.electric_count))
        self.electric=(hstack((base.er,zero_scalar),format='csr'),
            hstack((zero_meridional,toroidal),format='csr'),hstack((base.ez,zero_scalar),format='csr'))
        self.curl=(hstack((zero_meridional,curl_r),format='csr'),
            hstack((base.curl,zero_scalar),format='csr'),hstack((zero_meridional,curl_z),format='csr'))
        nr,dr,_,_=bases(.5*mesh.rf**2,np.array([.5*mesh.radius**2]),degree)
        nv,dv,_,_=bases(mesh.zf,mesh.zf[[0,-1]],degree)
        self.poloidal_null=kron(trace_kernel(np.vstack((nr.toarray(),dr.toarray()))),
            trace_kernel(np.vstack((nv.toarray(),dv.toarray()))),format='csr')
        self.toroidal_magnetic_count=base.magnetic_count;npoloidal=self.poloidal_null.shape[1]
        if npoloidal<1:raise ValueError('Enough independent nodes for compact poloidal value and derivative traces required')
        zero_toroidal=csr_matrix((points,base.magnetic_count));zero_poloidal=csr_matrix((points,npoloidal))
        self.magnetic=(hstack((zero_toroidal,curl_r@self.poloidal_null),format='csr'),
            hstack((base.magnetic,zero_poloidal),format='csr'),
            hstack((zero_toroidal,curl_z@self.poloidal_null),format='csr'))
        self.flux=hstack((zero_toroidal,scalar.multiply(self.reference_r[:,None]**2)@self.poloidal_null),format='csr')
        self.gradient=np.vstack((base.gradient,np.zeros((self.scalar_count,self.scalar_count))))
        self.electric_count=base.electric_count+self.scalar_count
        self.magnetic_count=base.magnetic_count+npoloidal
