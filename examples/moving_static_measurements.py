"""Independent full-space sphere energy, generalized force and vacuum controls."""
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_geometry import MovingDomain
from plasma_bench.moving_plasma import MovingPlasma,MovingState
from .moving_cases import provenance
from .verification_records import write_record


def sphere(n,scale=1.,shift=0.,matching=2.,modes=8,quadrature=5):
    d=MovingDomain(n,2*n,matching_radius=matching,vacuum_layers=round(n*(matching-1)),
                   modes=modes,quadrature=quadrature)
    b=np.sqrt(MU0);s=MovingPlasma(d,.5*b,pressure=0.,bias=.3*b)
    x=scale*d.fluid_reference;x[:,1]+=shift
    e=s.evaluate(MovingState(x,np.zeros_like(x)))
    exact=np.pi*(1/scale-.6*scale+.09*scale**3)
    derivative=np.pi*(-1/scale**2-.6+.27*scale**2)
    measured=-np.sum(e["force"]*d.fluid_reference)
    return {"radial_cells":n,"angular_cells":2*n,"radius_m":1.,"scale":scale,"axial_shift_m":shift,
            "matching_radius_m":matching,"vacuum_layers":round(n*(matching-1)),"modes":modes,
            "quadrature":quadrature,"initial_plasma_B_t":float(b),"bias_t":float(.3*b),
            "exact_energy_j":float(exact),"measured_energy_j":e["magnetic"].energy,
            "energy_relative_error":float(e["magnetic"].energy/exact-1),
            "exact_denergy_dscale_j":float(derivative),"measured_denergy_dscale_j":float(measured),
            "denergy_dscale_relative_error":float(measured/derivative-1),
            "spurious_axial_force_n":float(e["force"][:,1].sum()),
            "vacuum_relative_residual":e["magnetic"].vacuum_residual,
            "physical_minus_weighted_energy_j":e["magnetic"].physical_energy-e["magnetic"].energy}


def measure(directory="outputs/moving_verification"):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=provenance()
    rows={"provenance":source,"completed":True,
          "interpretation":"Electromagnetic virtual-force fixture; Bn and tangential-field jump do not define a free-gas equilibrium",
          "homologous":[sphere(n,scale) for scale in (.9,1.,1.1) for n in (4,8,16,32)],
          "translation":[sphere(n,shift=.3) for n in (4,8,16,32)],
          "modes":[sphere(12,shift=.4,modes=m) for m in (1,2,4,8,12)],
          "matching":[sphere(12,shift=.4,matching=r) for r in (1.5,2.,3.)],
          "quadrature":[sphere(12,shift=.4,quadrature=q) for q in (3,5,7)]}
    if provenance()!=source:
        raise RuntimeError("Source changed during static measurements")
    write_record(directory/"static.json",rows)
    print("Finished static moving-field studies",flush=True)
    return rows


if __name__=="__main__":
    measure()
