"""Independent clipped-polygon moments and discontinuous common-support integrals."""
import numpy as np
from examples.material_overlap import intersection,intersections,quadrature
from plasma_bench.moving_geometry import MovingDomain


def test_polygon_intersection_has_exact_area_and_radial_moment():
    triangle=np.array([[0.,0.],[2.,0.],[0.,2.]])
    box=np.array([[.5,.25],[1.5,.25],[1.5,1.25],[.5,1.25]])
    poly=intersection(triangle,box)
    values=np.zeros(2)
    for x,w in quadrature([poly],4):values+=np.array([w.sum(),np.sum(w*x[:,0])])
    # Integrate x from 1/2 to 3/4 with unit height, then height 7/4-x to 3/2.
    np.testing.assert_allclose(values,[23/32,83/128],atol=3e-15)


def test_overlap_integrates_exact_discontinuous_field_and_graded_volume():
    def record(polygons):
        points=np.asarray(polygons).reshape(-1,2)
        return {"t0_reference_coordinates":points,"t0_cells":np.arange(len(points)).reshape(-1,4)}
    a=record([[[.2,-.5],[1.,-.5],[1.,.5],[.2,.5]]])
    b=record([[[.2,-.5],[.83,-.5],[.83,.5],[.2,.5]],[[.83,-.5],[1.,-.5],[1.,.5],[.83,.5]]])
    mass=0.
    for x,w in quadrature(intersections(a,b),3):mass+=np.sum(2*np.pi*x[:,0]*w*(x[:,0]>.83))
    np.testing.assert_allclose(mass,np.pi*(1-.83**2),atol=4e-15)
    d=MovingDomain(4,8,radial_grading=1.5);f=MovingDomain(7,16,radial_grading=2.)
    a={"t0_reference_coordinates":d.fluid_reference,"t0_cells":d.cells[:d.fluid_cell_count]}
    b={"t0_reference_coordinates":f.fluid_reference,"t0_cells":f.cells[:f.fluid_cell_count]}
    mass=sum(float(np.sum(2*np.pi*x[:,0]*w)) for x,w in quadrature(intersections(a,b),4))
    edge=d.fluid_reference[d.boundary];r0=edge[:-1,0];r1=edge[1:,0];dz=np.diff(edge[:,1])
    exact=np.sum(np.pi/3*dz*(r0*r0+r0*r1+r1*r1))
    np.testing.assert_allclose(mass,exact,rtol=3e-14)
