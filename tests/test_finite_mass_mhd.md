# Complete coupled derivative and step checks

Original moving parent-species manufactured fields check all five conservative
rates and both pressure rates across three meshes. The changing-density current
energy is included in inferred pi_t; its substantial coarse error is required
to decrease, not hidden by a total-energy assertion. The same actual mass flux
must commute to dual n_t including numerical diffusion.

Uniform translation with axial bias checks stationary background preservation.
An unforced nonuniform full-current state conserves mass, axial/angular momenta
and total energy across actual RK steps. A deliberately excessive energy sink
must reject before any accepted field, time or budget changes. These are bounded
verification tests; no global stability, shock partition or formation claim.
