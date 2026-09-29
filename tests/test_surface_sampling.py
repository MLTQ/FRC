"""Independent affine/polynomial reconstruction across distinct graded meshes."""
import numpy as np
import pytest
from plasma_bench.moving_geometry import MovingDomain
from examples.surface_sampling import fields,nodal


@pytest.mark.parametrize("n,nt,grading",[(5,12,1.5),(8,20,2.)])
def test_common_physical_reference_points_reconstruct_analytic_fields(n,nt,grading):
    d=MovingDomain(n,nt,radial_grading=grading);x0=d.fluid_reference;x=x0.copy()
    x[:,0]*=.9;x[:,1]=1.1*x[:,1]+.05
    velocity=np.column_stack((.3*x0[:,0],-.1*x0[:,0]+.2*x0[:,1]))
    a=1+.2*x0[:,0]-.1*x0[:,1];q=2+3*x0[:,0]+4*x0[:,1]
    nc=d.fluid_cell_count
    arrays={"t0_reference_coordinates":x0,"t0_cells":d.cells[:nc],"t_coordinates":x,
            "t_velocity":velocity,"t_total_material_a":a,"t_current_over_radius":q,
            "t_density":np.ones(nc),"t_pressure":2*np.ones(nc),"t_volume":3*np.ones(nc),
            "t_cell_ohmic_heat":12*np.ones(nc)}
    config={"radial_cells":n,"angular_cells":nt}
    rng=np.random.default_rng(4);points=np.column_stack((rng.uniform(.02,.7,200),rng.uniform(-.5,.5,200)))
    values=fields(arrays,config,"t",points);r0,z0=points.T;r=.9*r0
    exact_position=np.column_stack((r,1.1*z0+.05))
    exact_b=np.column_stack((r*.1/1.1,2*(1+.2*r0-.1*z0)+r*.2/.9))
    np.testing.assert_allclose(values["position"],exact_position,atol=3e-15)
    np.testing.assert_allclose(values["magnetic"],exact_b,atol=2e-13)
    np.testing.assert_allclose(values["current"],r*(2+3*r0+4*z0),atol=2e-14)
    np.testing.assert_allclose(values["velocity"],np.column_stack((.3*r0,-.1*r0+.2*z0)),atol=3e-15)
    np.testing.assert_allclose(values["volume_ratio"],.9**2*1.1,atol=2e-14)
    np.testing.assert_array_equal(values["heat_density"],4.)
    np.testing.assert_allclose(nodal(arrays,config,"t",x0,"coordinates"),x,atol=3e-15)
    with pytest.raises(ValueError,match="outside"):
        nodal(arrays,config,"t",np.array([[1.1,0.]]),"velocity")
