"""Reproduce the extra time control requested by the completed fine-mesh evidence."""
from .resistive_records import ResistiveCase
from .resistive_dynamic_measurements import measure
from .surface_measurements import provenance


def case():
    return ResistiveCase("source_g2_32_cfl_quarter","source",32,radial_grading=2.,
                         cfl=.0375,max_timestep=.000125)


if __name__=="__main__":
    measure(case(),"outputs/surface_verification/dynamic",fingerprint=provenance)
