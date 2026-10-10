"""Chunk equivalence and complete discrete-wall identity controls."""
import numpy as np
from examples.strong_wall_reference import initialize
from examples.current_transport_sampling import TransportSampler
from examples.radial_term_sampling import radial_slice,velocity_rate,negative_curl,raw_faraday,projected_faraday,wall_projection
from plasma_bench.kernel_induction import difference_curl
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates


def test_chunked_original_fields_and_complete_wall_projection():
    _,d,u,s,b=initialize(dict(kind='contact',elements=2,points=16,field_nr=8,field_nz=16))
    sample=TransportSampler(d,s.geometry.mesh,8);part=radial_slice(sample,11,31)
    c=ModalWallCoordinates(s,u,b);du,db=c.unpack(1e-5*np.random.default_rng(392).normal(size=c.size))
    state=u+du;rate=du*5
    assert np.allclose(velocity_rate(part,state,rate),sample.fluid(state,rate)[4][11:31],rtol=1e-13,atol=1e-13)
    electric=np.random.default_rng(920).normal(size=s.geometry.faraday.shape[1])
    assert np.allclose(part.electric(electric),sample.electric(electric)[11:31],rtol=1e-13,atol=1e-13)
    assert np.array_equal(raw_faraday(s.geometry.mesh,electric),difference_curl(s.geometry,electric))
    p=wall_projection(s.geometry.mesh)
    assert np.array_equal(p,s.geometry.radial_projection)
    assert np.array_equal(projected_faraday(s.geometry.mesh,electric,p),s.geometry.project_curl(difference_curl(s.geometry,electric)))
    assert np.allclose(negative_curl(part,electric),part.magnetic(raw_faraday(s.geometry.mesh,electric)),rtol=2e-11,atol=2e-9)
