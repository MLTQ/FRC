# Sampling invariance controls

Checks actual conservative velocity differentiation and all electric components
in a radial chunk against unsliced independent physical reconstruction. Requires
bitwise identity of full raw incidence and wall-projected coefficients with
the original operator, including exact equality of auxiliary-axial and full
radial projection matrices. The independently differentiated physical strong
curl must agree with the reconstructed unprojected Faraday representation.
No accuracy improvement is inferred from these identity checks; analytic
parent-term checks and full physical refinements supply separate evidence.
