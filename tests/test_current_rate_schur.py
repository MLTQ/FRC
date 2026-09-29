"""Compare algebraic elimination to independent original mixed direct solves."""
import numpy as np
import pytest
from examples.current_rate_measurements import discrete as mixed
from examples.current_rate_schur import discrete as schur
from examples.current_thermal_cases import protocol
from examples.current_thermal_evolution import initialize


@pytest.mark.parametrize('nr', [4,8,16])
def test_positive_rate_elimination_preserves_both_variable_density_forcings(nr):
    case=next(c for c in protocol() if c['name']=='original_8');case.update(nr=nr,nz=2*nr)
    _,model,psi,g,ue,ui,_=initialize(case)
    direct,dr=mixed(model,psi,g,ue,ui);reduced,sr=schur(model,psi,g,ue,ui)
    assert max(dr+sr)<1e-10
    for a,b in zip(direct,reduced):
        assert np.linalg.norm(a-b)<1e-10*np.linalg.norm(a)
