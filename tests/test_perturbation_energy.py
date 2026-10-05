"""Independent high-precision conservative-EOS curvature verification."""
from decimal import Decimal, localcontext
import numpy as np
import pytest
from examples.perturbation_energy_reference import thermal_quadratic


@pytest.mark.parametrize('gamma,a,b',[('1.6666666666666666666666666666666667','.27','-.43'),('1.4','-.31','.19'),('1.01','.7','-.007')])
def test_thermal_hessian_matches_original_conservative_eos(gamma,a,b):
    with localcontext() as context:
        context.prec=70
        g=Decimal(gamma);a=Decimal(a);b=Decimal(b);s=Decimal('.37');h=Decimal('1e-10')
        def original(sign):
            n=1+sign*h*a;C=s+sign*h*(b+s*a)
            return (g*n.ln()+C/n).exp()
        expected=(original(1)+original(-1)-2*s.exp())/(2*h*h*s.exp())
    actual=thermal_quadratic(1.,float(a),float(b),float(g))
    assert np.isclose(actual,float(expected),rtol=2e-14,atol=0.)
    assert actual>0


def test_quadratic_energy_uses_conservative_density_entropy_cross_work():
    # At fixed entropy density C and s=2, delta s=-2*delta n/n.
    # The original EOS U=n^(3/2)*exp(2/n) has U''(1)/2 = 11*exp(2)/8.
    a=.125
    value=thermal_quadratic(np.exp(2.),a,-2*a,1.5)
    assert np.isclose(value,11*np.exp(2.)*a*a/8,rtol=2e-15)
    with pytest.raises(ValueError):thermal_quadratic(1.,a,0.,1.)
