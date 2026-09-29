"""Independent interpolation and conservative-restriction checks for study norms."""
import numpy as np
from plasma_bench.moving_geometry import MovingDomain
from examples.resistive_analysis import material_samples,restrict_cells


def test_actual_field_sampling_reproduces_linear_a_on_nested_deformed_sectors():
    n=4;nf=8;d=MovingDomain(nf,2*nf);x=d.fluid_reference.copy()
    x[:,0]=1.1*d.fluid_reference[:,0]
    x[:,1]=.9*d.fluid_reference[:,1]+.02*d.fluid_reference[:,0]+.05
    c,alpha,beta=.001,.0002,-.0003
    arrays={"t0_coordinates":x,"t0_cells":d.cells[:d.fluid_cell_count],
            "t0_total_material_a":c+alpha*x[:,0]+beta*x[:,1],
            "t0_current_over_radius":.3+.2*x[:,0]-.1*x[:,1]}
    magnetic,current,weight=material_samples(arrays,"t0",n,nf,subdivisions=2)
    g,_=np.polynomial.legendre.leggauss(5);nodes=np.r_[(g+1)/4,.5+(g+1)/4]
    xi,tau=np.meshgrid(nodes,nodes,indexing="ij");xi=xi.ravel();tau=tau.ravel()
    ir,it=np.divmod(np.arange(2*n*n),2*n)
    shell=(ir[:,None]+xi[None,:])/n
    j=2*it[:,None]+np.floor(2*tau)[None,:];fraction=(2*tau)%1
    angle0=j*np.pi/(2*nf);angle1=(j+1)*np.pi/(2*nf)
    r=shell*((1-fraction)*np.sin(angle0)+fraction*np.sin(angle1))
    z=-shell*((1-fraction)*np.cos(angle0)+fraction*np.cos(angle1))
    rr=1.1*r;zz=.9*z+.02*r+.05
    expected=np.stack((-rr*beta,2*c+3*alpha*rr+2*beta*zz),axis=-1)
    np.testing.assert_allclose(magnetic,expected,rtol=2e-13,atol=1e-17)
    np.testing.assert_allclose(current,rr*(.3+.2*rr-.1*zz),rtol=2e-14,atol=1e-15)
    boundary=x[d.boundary];r0,r1=boundary[:-1,0],boundary[1:,0]
    exact_volume=np.pi/3*np.sum(np.diff(boundary[:,1])*(r0*r0+r0*r1+r1*r1))
    assert abs(weight.sum()/exact_volume-1)<2e-14


def test_cell_restriction_preserves_nonuniform_integrated_heat_and_pressure():
    n,nf=2,4;volume=np.arange(1,2*nf*nf+1,dtype=float)
    pressure=np.arange(len(volume),dtype=float)**2;heat=np.sqrt(volume)
    arrays={"t0_volume":volume,"t0_pressure":pressure,"t0_cell_ohmic_heat":heat}
    for key,values,integrated in (("pressure",pressure,False),("cell_ohmic_heat",heat,True)):
        actual=restrict_cells(arrays,"t0",n,nf,key);expected=[]
        for i in range(n):
            for j in range(2*n):
                ids=[(2*i+di)*2*nf+2*j+dj for di in (0,1) for dj in (0,1)]
                expected.append(sum(values[k]*(1 if integrated else volume[k]) for k in ids)/sum(volume[k] for k in ids))
        np.testing.assert_allclose(actual,expected,rtol=2e-15)
