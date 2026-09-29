# Audited physical thermal report

Requires all18 rate and9 transport cases with matching source and audit/raw hashes.
Six panels distinguish matched old/new nodal axis rates, continuous physical
proxy rates, pressure accuracy, full-transit shear error, isothermal false
battery and conservation. Relative/absolute and nodal/physical norms stay labeled.
No plot or table converts a conservation result into a general accuracy claim.

`full_field_conditioning` adds physical magnetic-energy-norm and T/s RMS errors
for Faraday replay of the rounded full E array, separate from the stiff solve's
reduced representation. It uses the already verified physical magnetic Hessian;
this is a conditioning diagnostic, not an independent continuum oracle. Exact-
zero rates have no fabricated relative denominator. `time_control` independently
samples the two final time-step-control states against the same entropy signal.

The report helper and companion have their own archive/hashes outside the frozen
calculation closure. The figure needs visual inspection. All raw controls and
signed energy/proxy drift stay in the report, including limitations. The result
does not establish full moving plasma, physical transport or formation.
