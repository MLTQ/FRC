"""Physical volume moments and independent moving-gas accuracy on graded meshes."""
import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from plasma_bench.moving_geometry import MovingDomain
from examples.surface_references import expansion_case


@pytest.mark.parametrize("grading",[1.5,2.])
def test_graded_volume_and_second_moments_match_revolved_polygon(grading):
    d=MovingDomain(7,18,radial_grading=grading)
    e=d.elements(d.fluid_reference)[0];n=d.fluid_cell_count
    r,z=e.points[:n,:,0],e.points[:n,:,1];w=2*np.pi*r*e.weights[:n]
    measured=np.array([w.sum(),np.sum(w*r*r),np.sum(w*z*z)])
    # Independent disks along z, with a linear radius on each polygon edge.
    edge=d.fluid_reference[d.boundary];q,weights=leggauss(8);exact=np.zeros(3)
    for left,right in zip(edge[:-1],edge[1:]):
        z=(left[1]+right[1])/2+(right[1]-left[1])*q/2
        r=(left[0]+right[0])/2+(right[0]-left[0])*q/2
        dz=weights*(right[1]-left[1])/2
        exact+=np.array([np.sum(dz*np.pi*r*r),np.sum(dz*np.pi*r**4/2),np.sum(dz*np.pi*r*r*z*z)])
    np.testing.assert_allclose(measured,exact,rtol=3e-14,atol=1e-14)


@pytest.mark.parametrize("grading",[1.5,2.])
def test_graded_actual_gas_expansion_refines_against_continuum(grading):
    values=[]
    for n in (4,8,16):
        row,_=expansion_case(n,grading,end_time=.02)
        values.append(row["errors"]["velocity_relative_mass_L2"])
        assert row["final"]["mass_kg"]==row["initial"]["mass_kg"]
        assert abs(row["final"]["response_energy_residual_j"])<2e-7
    assert np.all(np.asarray(values[:-1])/values[1:]>2)
    assert values[-1]<.04


def test_invalid_or_collapsed_radial_mesh_is_rejected():
    for grading in (0.,-1.,np.nan,np.inf,1e4):
        with pytest.raises(ValueError,match="grading"):
            MovingDomain(8,16,radial_grading=grading)
