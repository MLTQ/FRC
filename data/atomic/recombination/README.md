# Ground-state argon radiative recombination inputs

Raw sources are retained unchanged; `argon_rr.json` stores selected rows, their
source line numbers, physical units, URLs and SHA256 hashes. The complete
ladder covers initial charges 1..18. Binding energies remain in the separate
NIST-based `../argon/argon.json`; fitting parameters do not replace them.

- `badnell_clist_ev.dat`: [author's total RR table](https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/clist_eV),
  header dated 20230511. Select Z=18, M=1 and initial charge q=Z-N.
  It supplies q=3..18. T0,T1,T2 are eV and A is cm³/s.
  Use the [published fit and data conventions](https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/).
  Nominal range is q²(10..10^7) K. Initial ground-state rates only are used.
  [Badnell 2006](https://arxiv.org/abs/astro-ph/0604144) describes bare through
  Na-like ions; [sequence-specific publications](https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/pub.html)
  cover subsequent extensions. These are theoretical rates, not experimental
  validation of this plasma model.
- `verner_pl.dat` and `verner_pl.txt`: [archival power-law table](https://www.pa.uky.edu/~verner/rec.html).
  Select Z=18 and recombined electron counts N=18,17, corresponding to initial
  charges q=Z-N+1=1,2. Their coefficients are respectively A=3.77e-13,
  1.95e-12 cm³/s and exponents .651,.752 in A*(T/10000 K)^(-exponent).
  Verner's [routine attribution](https://www.pa.uky.edu/~verner/dima/rec//rrfit.f)
  identifies these argon rates with Shull & Van Steenberg (1982).
  The distributed table does not specify a validity interval. It is stored as
  **unknown**, not invented or silently treated as valid. These legacy fits
  dominate the RR events in our low-charge source examples; improved data is
  still required for a credible preparation model.
- `badnell_cfout_ev.gz`: [original author calculation and fit comparison](https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/cfout_eV.gz).
  `argon_rr_original.json` extracts its Z=18,M=1 rows with original columns
  unchanged: temperature, original rate, fitted rate, absolute relative error.
  The downloaded master file includes 12 argon ground-state stages; q=8,16,17,18
  are absent. Original-rate comparisons cover those 12 stages only. Values and
  temperatures in this file are rounded, unlike the underlying calculation.

These data describe zero-density RR, assuming Maxwellian electrons. They do
not include dielectronic/three-body recombination, density suppression,
metastable populations or photon reabsorption. The model's photon-energy
moment differentiates the rate fits, an additional approximation that their
rate-value accuracy does not by itself validate. Immediate cascade to ground
and unrestricted photon escape are explicit source-study assumptions.
No restricted OPEN-ADAS files are used. These are local research archives;
no remote publication is configured or requested.
